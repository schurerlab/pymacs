import base64
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from pymacs_cgenff_funnel import parameterize_mol2


class CGenFFFunnelClientTests(unittest.TestCase):
    def test_parameterize_posts_mol2_and_decodes_stream(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            token = folder / "token"
            mol2 = folder / "DR7.cgenff.mol2"
            token.write_text("x" * 32, encoding="utf-8")
            mol2.write_text("@<TRIPOS>MOLECULE\nDR7\n", encoding="utf-8")
            response = MagicMock()
            response.read.return_value = json.dumps({
                "str_b64": base64.b64encode(b"* CGenFF stream\n").decode("ascii"),
                "log_b64": base64.b64encode(b"done\n").decode("ascii"),
            }).encode("utf-8")
            response.__enter__.return_value = response
            with patch("pymacs_cgenff_funnel.urlopen", return_value=response) as mocked:
                stream, log = parameterize_mol2(
                    endpoint="https://kyle.example.ts.net/", token_file=token,
                    ligand_code="dr7", mol2_path=mol2,
                )
            self.assertEqual(stream, "* CGenFF stream\n")
            self.assertEqual(log, "done\n")
            request = mocked.call_args.args[0]
            self.assertEqual(request.full_url, "https://kyle.example.ts.net/v1/parameterize")
            self.assertEqual(request.get_header("Authorization"), "Bearer " + "x" * 32)
