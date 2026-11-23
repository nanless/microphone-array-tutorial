"""Independent byte and image controls for the publication review repairs."""

import hashlib
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch
import zlib

from PIL import Image, PngImagePlugin

from scripts import quality_check


class FigureIdentityTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.script = self.root / 'draw.py'
        self.script.write_text('# independently authored fixture\n')
        self.path = self.root / 'fig01_example.png'
        metadata = PngImagePlugin.PngInfo()
        metadata.add_text('SourceScript', 'draw.py')
        metadata.add_text('SourceScriptDigest', hashlib.sha256(self.script.read_bytes()).hexdigest())
        metadata.add_text('FigureFilename', self.path.name)
        Image.new('RGB', (16, 16), (20, 40, 60)).save(self.path, pnginfo=metadata)
        context = patch.object(quality_check, 'ROOT', self.root)
        context.start()
        self.addCleanup(context.stop)

    def test_valid_decoded_png_passes_without_changing_its_bytes(self):
        before = self.path.read_bytes()
        self.assertEqual(quality_check.png_provenance_issues(self.path, self.script), [])
        self.assertEqual(self.path.read_bytes(), before)

    def test_another_valid_png_from_the_same_producer_cannot_supply_this_figure(self):
        other = self.root / 'fig02_other.png'
        other.write_bytes(self.path.read_bytes())
        issues = quality_check.png_provenance_issues(other, self.script)
        self.assertTrue(any('FigureFilename' in issue for issue in issues))

    def test_valid_chunk_crc_does_not_hide_undecodable_pixel_data(self):
        raw = bytearray(self.path.read_bytes())
        start = raw.index(b'IDAT')
        length = struct.unpack_from('>I', raw, start - 4)[0]
        raw[start + 4] = 0  # invalid deflate header, with a matching outer PNG CRC
        checksum = zlib.crc32(raw[start:start + 4 + length])
        struct.pack_into('>I', raw, start + 4 + length, checksum)
        self.path.write_bytes(raw)
        with Image.open(self.path) as image:
            image.verify()  # outer chunk structure alone still succeeds
        with self.assertRaises(OSError):
            quality_check.png_provenance_issues(self.path, self.script)


if __name__ == '__main__':
    unittest.main()
