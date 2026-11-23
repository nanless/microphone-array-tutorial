"""Bind a printed PDF to its actual dated HTML input.

These are finite checks around a local build. They do not make concurrent
filesystem access atomic or detect every change-and-restore sequence.
"""
from __future__ import annotations

from dataclasses import dataclass
import datetime
import hashlib
from html.parser import HTMLParser
from pathlib import Path
import re

from codes.chapters.ch00.io_contracts import validate_parent_chain


@dataclass(frozen=True)
class PDFInputIdentity:
    source_digest: str
    html_sha256: str
    build_date: str

    def __post_init__(self):
        if not isinstance(self.source_digest, str) or not re.fullmatch(r'[0-9a-f]{12}', self.source_digest):
            raise ValueError('HTML source digest must be twelve lowercase hexadecimal characters')
        if not isinstance(self.html_sha256, str) or not re.fullmatch(r'[0-9a-f]{64}', self.html_sha256):
            raise ValueError('HTML input SHA256 must be sixty-four lowercase hexadecimal characters')
        if not isinstance(self.build_date, str) or not re.fullmatch(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', self.build_date):
            raise ValueError('HTML build date must be YYYY-MM-DD')
        datetime.date.fromisoformat(self.build_date)

    def metadata(self):
        return {'/SourceDigest': self.source_digest,
                '/HTMLInputSHA256': self.html_sha256,
                '/BuildDate': self.build_date}

    @classmethod
    def from_metadata(cls, metadata):
        return cls(metadata.get('/SourceDigest'), metadata.get('/HTMLInputSHA256'),
                   metadata.get('/BuildDate'))


class _HeadIdentityParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.in_head = False
        self.heads = 0
        self.closed_heads = 0
        self.values = {'tutorial-source-digest': [], 'tutorial-build-date': []}

    def handle_starttag(self, tag, attrs):
        if tag == 'head':
            self.heads += 1
            self.in_head = True
        if tag != 'meta' or not self.in_head:
            return
        names = [value for key, value in attrs if key == 'name']
        if not any(name in self.values for name in names):
            return
        attributes = dict(attrs)
        if len(attributes) != len(attrs) or len(names) != 1:
            raise ValueError('duplicate HTML provenance attributes')
        self.values[names[0]].append(attributes.get('content'))

    def handle_endtag(self, tag):
        if tag == 'head':
            if not self.in_head:
                raise ValueError('unbalanced HTML head')
            self.closed_heads += 1
            self.in_head = False


def html_input_identity(html_bytes: bytes) -> PDFInputIdentity:
    if not isinstance(html_bytes, bytes):
        raise TypeError('HTML input must be captured as bytes')
    parser = _HeadIdentityParser()
    parser.feed(html_bytes.decode('utf-8'))
    parser.close()
    if parser.heads != 1 or parser.closed_heads != 1 or parser.in_head:
        raise ValueError('HTML input needs one complete head')
    if any(len(values) != 1 for values in parser.values.values()):
        raise ValueError('HTML input needs unique source and build-date metadata')
    return PDFInputIdentity(parser.values['tutorial-source-digest'][0],
                            hashlib.sha256(html_bytes).hexdigest(),
                            parser.values['tutorial-build-date'][0])


def _canonical_html(builder, build_date):
    result = builder(build_date)
    text = result[0] if isinstance(result, tuple) else result
    if not isinstance(text, str):
        raise TypeError('canonical HTML builder must return text')
    return text.encode('utf-8')


def verify_canonical_html_input(html_bytes, builder, source_digest) -> PDFInputIdentity:
    """Check captured bytes against current sources using their explicit date."""
    before = source_digest()
    identity = html_input_identity(html_bytes)
    if identity.source_digest != before:
        raise ValueError('HTML input source digest is stale')
    canonical = _canonical_html(builder, identity.build_date)
    if canonical != html_bytes:
        raise ValueError('HTML input differs from the canonical body for its build date')
    if source_digest() != before:
        raise ValueError('publication sources changed during HTML validation')
    return identity


def require_unchanged_pdf_input(identity, html_path, source_digest):
    path = validate_parent_chain(html_path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != identity.html_sha256:
        raise ValueError('private PDF HTML input changed after capture')
    if source_digest() != identity.source_digest:
        raise ValueError('publication sources changed after PDF HTML capture')


def verify_pdf_input_metadata(reader, builder, source_digest) -> PDFInputIdentity:
    """Rebuild dated HTML in memory, independently of today's combined file."""
    identity = PDFInputIdentity.from_metadata(reader.metadata or {})
    before = source_digest()
    if identity.source_digest != before:
        raise ValueError('PDF input source digest is stale')
    canonical = _canonical_html(builder, identity.build_date)
    if html_input_identity(canonical) != identity:
        raise ValueError('PDF metadata differs from its canonical dated HTML input')
    if source_digest() != before:
        raise ValueError('publication sources changed during PDF input validation')
    return identity
