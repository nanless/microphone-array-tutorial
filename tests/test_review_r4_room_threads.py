"""Room thread execution and metadata contracts; no published asset writes."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from codes.chapters.appendix_b.examples import check_room_assets as audit
from codes.chapters.appendix_b.examples import room_srp_exercise as room


class _Constants:
    def __init__(self, threads):
        self.threads = threads

    def get(self, name):
        assert name == 'num_threads'
        return self.threads

    def set(self, name, value):
        assert name == 'num_threads'
        self.threads = value


class RoomThreadContractsTest(unittest.TestCase):
    def test_private_room_factory_sets_actual_threads_without_changing_env(self):
        observed = []
        constants = _Constants(10)

        class ObservedRoom:
            def __init__(self, *args, **kwargs):
                pass

            def set_sound_speed(self, value):
                pass

            def add_source(self, value):
                pass

            def add_microphone_array(self, value):
                pass

            def compute_rir(self):
                observed.append(constants.get('num_threads'))

        pra = SimpleNamespace(constants=constants, ShoeBox=ObservedRoom,
                              Material=lambda absorption: absorption)
        for requested in (None, '10'):
            with self.subTest(requested=requested), patch.dict(os.environ):
                if requested is None:
                    os.environ.pop('PRA_NUM_THREADS', None)
                else:
                    os.environ['PRA_NUM_THREADS'] = requested
                before = dict(os.environ)
                constants.threads = 10
                room._make_room(pra, [6., 6., 3.], [[6., 5., 3.]] * 4, .4, 40)
                self.assertEqual(observed[-1], 2)
                self.assertEqual(dict(os.environ), before)

    def test_runtime_records_requested_and_effective_threads_separately(self):
        # Keep this metadata unit test usable in the NumPy-only environment.
        with tempfile.TemporaryDirectory() as temporary, patch.dict(
                'sys.modules', {'scipy': SimpleNamespace(__version__='fixture')}):
            library = Path(temporary)
            for name in ('room.py', 'acoustics.py', 'parameters.py'):
                (library / name).write_text('# runtime identity fixture\n', encoding='utf-8')
            pra = SimpleNamespace(__file__=str(library / '__init__.py'),
                                  __version__=room.PRA_VERSION, constants=_Constants(2))
            for requested in (None, '10'):
                with self.subTest(requested=requested), patch.dict(os.environ):
                    if requested is None:
                        os.environ.pop('PRA_NUM_THREADS', None)
                    else:
                        os.environ['PRA_NUM_THREADS'] = requested
                    environment = room.runtime_environment(pra)
                    self.assertEqual(environment['pra_num_threads'], requested)
                    self.assertIs(type(environment['pra_num_threads_effective']), int)
                    self.assertEqual(environment['pra_num_threads_effective'], 2)
                    self.assertEqual(os.environ.get('PRA_NUM_THREADS'), requested)

    def test_metadata_checks_effective_type_and_cross_document_binding(self):
        # This is a metadata fixture, not an assertion of a fresh PRA replay.
        # Actual original-library/WAV replay is performed separately.
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / 'room'
            shutil.copytree(audit.ROOM, directory)
            manifest = json.loads((directory / 'MANIFEST.json').read_text(encoding='utf-8'))
            result = json.loads((directory / 'RESULTS.json').read_text(encoding='utf-8'))
            for document in (manifest, result):
                document['source_sha256'] = room.source_hashes()
                document['environment']['pra_num_threads'] = None
                document['environment']['pra_num_threads_effective'] = 2
            result['generator']['sha256'] = room.source_hashes()[room.SOURCE_PATHS[0]]

            def write_fixture(m, r):
                payload = (json.dumps(m, ensure_ascii=False, indent=2) + '\n').encode('utf-8')
                (directory / 'MANIFEST.json').write_bytes(payload)
                r['assets']['audio_manifest']['sha256'] = hashlib.sha256(payload).hexdigest()
                (directory / 'RESULTS.json').write_text(
                    json.dumps(r, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

            write_fixture(copy.deepcopy(manifest), copy.deepcopy(result))
            self.assertFalse(audit.check_assets(directory)['validation']['rir_and_output_replayed'])
            for invalid in (None, True, 2., '2', 10, 'missing'):
                with self.subTest(effective=invalid):
                    m, r = copy.deepcopy(manifest), copy.deepcopy(result)
                    for document in (m, r):
                        if invalid == 'missing':
                            document['environment'].pop('pra_num_threads_effective')
                        else:
                            document['environment']['pra_num_threads_effective'] = invalid
                    write_fixture(m, r)
                    with self.assertRaisesRegex(ValueError, 'effective PRA threads'):
                        audit.check_assets(directory)
            for invalid in (2, False, ['2'], 'missing'):
                with self.subTest(raw_env=invalid):
                    m, r = copy.deepcopy(manifest), copy.deepcopy(result)
                    for document in (m, r):
                        if invalid == 'missing':
                            document['environment'].pop('pra_num_threads')
                        else:
                            document['environment']['pra_num_threads'] = invalid
                    write_fixture(m, r)
                    with self.assertRaisesRegex(ValueError, 'original PRA_NUM_THREADS'):
                        audit.check_assets(directory)
            m, r = copy.deepcopy(manifest), copy.deepcopy(result)
            m['environment']['pra_num_threads'] = '10'
            write_fixture(m, r)
            with self.assertRaisesRegex(ValueError, 'room environment binding'):
                audit.check_assets(directory)
            r['environment']['pra_num_threads'] = '10'
            write_fixture(m, r)
            self.assertFalse(audit.check_assets(directory)['validation']['rir_and_output_replayed'])


if __name__ == '__main__':
    unittest.main()
