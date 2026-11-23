"""Independent output and executed-source controls for chapter-owned audit CLIs."""
import contextlib
import hashlib
import importlib
import os
import py_compile
import sys
from types import SimpleNamespace
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import wave

import numpy as np

from codes.chapters.ch01.examples import audit_libmysofa_loudness as loudness
from codes.chapters.ch06 import aec_doubletalk_experiment as doubletalk


class AuditOutputBoundaryTests(unittest.TestCase):
    def test_loudness_rejects_formal_or_historical_report_before_compilation(self):
        for target in (loudness.LOCK, loudness.ROOT/'codes/chapters/ch01/reports/libmysofa_loudness.json'):
            with self.subTest(target=target), patch('sys.argv', ['audit', '--report', str(target)]), \
                 patch.object(loudness, 'run_audit') as run, self.assertRaises(ValueError):
                loudness.main()
            run.assert_not_called()

    def test_loudness_safe_external_report_and_symlink_parent(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); target = root/'report.json'
            with patch('sys.argv', ['audit', '--report', str(target)]), \
                 patch.object(loudness, 'run_audit', return_value={'expected_behavior_observed': True}), \
                 contextlib.redirect_stdout(io.StringIO()):
                loudness.main()
            self.assertTrue(target.is_file())
            alias=root/'alias';alias.symlink_to(root,target_is_directory=True)
            with patch('sys.argv', ['audit', '--report', str(alias/'bad.json')]), \
                 patch.object(loudness, 'run_audit') as run, self.assertRaises(ValueError):
                loudness.main()
            run.assert_not_called()

    def test_doubletalk_rejects_unsafe_output_before_original_run(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name); pair=root/'inputs';pair.mkdir();library=root/'native.dylib';library.write_bytes(b'fixture')
            existing=root/'existing.wav';existing.write_bytes(b'preserve')
            alias=root/'alias';alias.symlink_to(root,target_is_directory=True)
            targets=[doubletalk.ROOT/'new.wav', pair/'new.wav', existing, alias/'new.wav']
            for target in targets:
                argv=['audit','--speex-library',str(library),'--pair-dir',str(pair),'--output-wav',str(target)]
                with self.subTest(target=target), patch('sys.argv',argv), patch.object(doubletalk,'run') as run, \
                     self.assertRaises((ValueError,FileExistsError)):
                    doubletalk.main()
                run.assert_not_called()
            self.assertEqual(existing.read_bytes(),b'preserve')

    def test_doubletalk_new_external_wav_is_exclusive_and_preserves_samples(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);pair=root/'inputs';pair.mkdir();library=root/'native.dylib';library.write_bytes(b'fixture')
            target=root/'output'/'new.wav'
            argv=['audit','--speex-library',str(library),'--pair-dir',str(pair),'--output-wav',str(target)]
            with patch('sys.argv',argv), patch.object(doubletalk,'run',return_value=({},np.array([1,-2],dtype=np.int16))), \
                 contextlib.redirect_stdout(io.StringIO()):
                doubletalk.main()
            with wave.open(str(target),'rb') as wav:
                self.assertEqual((wav.getnchannels(),wav.getsampwidth(),wav.getframerate()),(1,2,16000))
                np.testing.assert_array_equal(np.frombuffer(wav.readframes(2),dtype='<i2'),[1,-2])
            with patch('sys.argv',argv), patch.object(doubletalk,'run') as run,self.assertRaises(FileExistsError):
                doubletalk.main()
            run.assert_not_called()

    def test_doubletalk_postflight_rejects_output_replaced_by_link(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);pair=root/'inputs';pair.mkdir();library=root/'native.dylib';library.write_bytes(b'fixture')
            victim=root/'victim';victim.write_bytes(b'preserve');target=root/'new.wav'
            def original(*args):
                target.symlink_to(victim)
                return {},np.array([1],dtype=np.int16)
            argv=['audit','--speex-library',str(library),'--pair-dir',str(pair),'--output-wav',str(target)]
            with patch('sys.argv',argv), patch.object(doubletalk,'run',side_effect=original),self.assertRaises(ValueError):
                doubletalk.main()
            self.assertEqual(victim.read_bytes(),b'preserve')


class LcmvNormalizationReviewTests(unittest.TestCase):
    def test_nonzero_response_cannot_disappear_during_column_normalization(self):
        from codes.chapters.ch05.core.beamforming import lcmv_weights
        for scale, response in ((1e308,1e-20),(1e200,1e-200)):
            with self.assertRaisesRegex(ValueError,'normalization loses'):
                lcmv_weights(np.eye(2),scale*np.eye(2),np.full(2,response))
        np.testing.assert_allclose(lcmv_weights(np.eye(3),1e200*np.eye(3),np.full(3,1e200)),np.ones(3))
        np.testing.assert_allclose(lcmv_weights(np.eye(2),np.eye(2),[1.,2.]),[1.,2.])
        np.testing.assert_array_equal(lcmv_weights(np.eye(2),np.eye(2),[0.,0.]),[0.,0.])


class OriginalPythonExecutionTests(unittest.TestCase):
    def test_verified_source_executes_instead_of_valid_old_pyc_without_changing_it(self):
        from codes.chapters.ch04.core import upstream_contracts as contracts
        package_name='review_original_bytecode_fixture'
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);package=root/package_name;package.mkdir();source=package/'__init__.py'
            source.write_text('value = 1\n');stamp=source.stat()
            cache=Path(py_compile.compile(
                str(source),doraise=True,
                invalidation_mode=py_compile.PycInvalidationMode.TIMESTAMP,
            ));old_cache=cache.read_bytes();old_mtime=cache.stat().st_mtime_ns
            source.write_text('value = 2\n');os.utime(source,ns=(stamp.st_atime_ns,stamp.st_mtime_ns))
            record={'python_files':{'source':{'path':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}}}
            sys.path.insert(0,str(root))
            try:
                normal=importlib.import_module(package_name);self.assertEqual(normal.value,1)
                sys.modules.pop(package_name)
                with contracts.source_imports(package_name,record):
                    actual=importlib.import_module(package_name)
                self.assertEqual(actual.value,2)
                self.assertEqual(actual.__micarray_original_source_sha256__,record['python_files']['source']['sha256'])
                self.assertEqual(cache.read_bytes(),old_cache);self.assertEqual(cache.stat().st_mtime_ns,old_mtime)
            finally:
                sys.path.remove(str(root));sys.modules.pop(package_name,None)

    def test_sourceless_bytecode_is_rejected_instead_of_executed(self):
        from codes.chapters.ch04.core import upstream_contracts as contracts
        package_name='review_sourceless_fixture'
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);package=root/package_name;package.mkdir();source=package/'__init__.py'
            source.write_text('value=7\n');cache=package/'__init__.pyc';py_compile.compile(str(source),cfile=str(cache),doraise=True);source.unlink()
            saved=cache.read_bytes();sys.path.insert(0,str(root))
            try:
                with contracts.source_imports(package_name),self.assertRaisesRegex(ValueError,'bytecode but no source'):
                    importlib.import_module(package_name)
                self.assertEqual(cache.read_bytes(),saved)
            finally:
                sys.path.remove(str(root));sys.modules.pop(package_name,None)

    def test_changed_source_is_rejected_before_its_body_runs(self):
        from codes.chapters.ch04.core import upstream_contracts as contracts
        package_name='review_rejected_source_fixture'
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);package=root/package_name;package.mkdir();source=package/'__init__.py';victim=root/'executed'
            source.write_text('from pathlib import Path\nPath('+repr(str(victim))+').write_text("bad")\n')
            record={'python_files':{'source':{'path':str(source),'sha256':'0'*64}}}
            sys.path.insert(0,str(root))
            try:
                with contracts.source_imports(package_name,record),self.assertRaisesRegex(ValueError,'execution bytes differ'):
                    importlib.import_module(package_name)
                self.assertFalse(victim.exists())
            finally:
                sys.path.remove(str(root));sys.modules.pop(package_name,None)

    def test_filterpy_same_filename_does_not_prove_preloaded_execution(self):
        from codes.chapters.ch09.examples import audit_tracking_upstream_interfaces as tracking
        with tempfile.TemporaryDirectory() as name:
            root=Path(name);package=root/'filterpy';package.mkdir();source=package/'__init__.py';source.write_text('value=1\n')
            identity={'checkout':str(root),'used_files':{'filterpy/__init__.py':{'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}}}
            prior={k:v for k,v in sys.modules.items() if k=='filterpy' or k.startswith('filterpy.')}
            for k in prior:sys.modules.pop(k)
            try:
                with patch.dict(sys.modules,{'filterpy':SimpleNamespace(__file__=str(source))}):
                    with self.assertRaisesRegex(ValueError,'lacks verified source-byte execution'):
                        tracking.loaded_filterpy_modules(identity)
            finally:
                for k in tuple(sys.modules):
                    if k=='filterpy' or k.startswith('filterpy.'):sys.modules.pop(k)
                sys.modules.update(prior)

    def test_said_executes_complete_verified_payload_without_loader_bytecode(self):
        from codes.chapters.ch04.examples import audit_said_compression as said
        if not (said.CACHE/'said/utils/compression.py').is_file():
            self.skipTest('optional pinned SAID source is unavailable')
        original_spec=said.importlib.util.spec_from_file_location
        loaders=[]
        def spec_without_executable_loader(name,path):
            result=original_spec(name,path)
            result.loader.exec_module=unittest.mock.Mock(side_effect=AssertionError('bytecode-capable loader used'))
            loaders.append(result.loader)
            return result
        with patch.object(said.importlib.util,'spec_from_file_location',side_effect=spec_without_executable_loader):
            result=said.build_report()
        for loader in loaders:loader.exec_module.assert_not_called()
        self.assertEqual(result['scope']['executed_source_sha256'],
                         hashlib.sha256((said.CACHE/'said/utils/compression.py').read_bytes()).hexdigest())

    def test_said_restores_temporary_module_entry_for_all_prior_states(self):
        from codes.chapters.ch04.examples import audit_said_compression as said
        if not (said.CACHE/'said/utils/compression.py').is_file():
            self.skipTest('optional pinned SAID source is unavailable')
        name='_masp_said_original_compression'
        absent=object()
        original=sys.modules.get(name,absent)
        try:
            for label,prior in [('None',None),('object',object()),('absent',absent)]:
                with self.subTest(prior_state=label):
                    if prior is absent:
                        sys.modules.pop(name,None)
                    else:
                        sys.modules[name]=prior
                    report=said.build_report()
                    self.assertIs(sys.modules.get(name,absent),prior)
                    self.assertEqual(name in sys.modules,prior is not absent)
                    self.assertEqual(report['scope']['executed_source_sha256'],
                                     hashlib.sha256((said.CACHE/'said/utils/compression.py').read_bytes()).hexdigest())
        finally:
            if original is absent:
                sys.modules.pop(name,None)
            else:
                sys.modules[name]=original

    def test_pra_preflight_failure_prevents_package_method_execution(self):
        from codes.chapters.ch04.examples import audit_upstream_doa as doa
        with patch.object(doa,'verify_checkout',return_value={}), \
             patch.object(doa.importlib.metadata,'version',return_value='0.10.0'), \
             patch.object(doa.contracts,'package_preflight',side_effect=ValueError('preflight sentinel')), \
             patch.object(doa,'_execute_report') as execute,self.assertRaisesRegex(ValueError,'preflight sentinel'):
            doa.build_report()
        execute.assert_not_called()


if __name__=='__main__':
    unittest.main()
