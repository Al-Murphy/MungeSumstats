import pathlib,subprocess,tempfile,unittest

SCRIPT=pathlib.Path(__file__).resolve().parents[2]/'scripts/benchmark-scale.py'
class VerifyIdentity(unittest.TestCase):
    def reject(self,mode,body):
        with tempfile.TemporaryDirectory() as name:
            d=pathlib.Path(name)
            panel=d/'panel.tsv';panel.write_text('chr\tpos\tref\talt\tid\taliases\n1\t1\tA\tC\trs1\t.\n1\t2\tG\tA\trs2\t.\n')
            source=d/'input.tsv';source.write_text('SNP\tCHR\tBP\tA1\tA2\tBETA\tFRQ\nrs1\t1\t1\tA\tC\t.125\t.25\nrs2\t1\t2\tG\tA\t.125\t.25\n')
            output=d/'output';output.mkdir()
            target=output/'normalized.tsv' if mode=='new' else output/'legacy.tsv'
            target.write_text(body)
            result=subprocess.run(['python3',str(SCRIPT),'verify','--panel',str(panel),'--input',str(source),'--rows','2','--mode',mode,'--output',str(output if mode=='new' else target),'--receipt',str(d/'receipt.json')],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0,result.stdout)
    def test_new_reordered_keys_do_not_match_source_rows(self):
        self.reject('new','source_row\tstatus\tchr\tpos\tref\talt\tid\teffect\teaf\n1\tmatch\t1\t2\tG\tA\trs2\t.125\t.25\n2\tmatch\t1\t1\tA\tC\trs1\t.125\t.25\n')
    def test_new_wrong_rsid_is_rejected(self):
        self.reject('new','source_row\tstatus\tchr\tpos\tref\talt\tid\teffect\teaf\n1\tmatch\t1\t1\tA\tC\trs2\t.125\t.25\n2\tmatch\t1\t2\tG\tA\trs2\t.125\t.25\n')
    def test_legacy_duplicate_identity_cannot_replace_another_source_row(self):
        self.reject('legacy','SNP\tCHR\tBP\tA1\tA2\tBETA\tFRQ\nrs1\t1\t1\tA\tC\t.125\t.25\nrs1\t1\t1\tA\tC\t.125\t.25\n')
    def test_typed_small_variant_cannot_be_relabeled_gene_burden(self):
        self.reject('new','source_row\tstatus\tchr\tpos\tref\talt\tid\teffect\teaf\tentity_type\tentity_id\n1\tmatch\t1\t1\tA\tC\trs1\t.125\t.25\tgene_burden\trs1\n2\tmatch\t1\t2\tG\tA\trs2\t.125\t.25\tsmall_variant\trs2\n')
    def test_typed_entity_id_must_match_its_original_source(self):
        self.reject('new','source_row\tstatus\tchr\tpos\tref\talt\tid\teffect\teaf\tentity_type\tentity_id\n1\tmatch\t1\t1\tA\tC\trs1\t.125\t.25\tsmall_variant\trs2\n2\tmatch\t1\t2\tG\tA\trs2\t.125\t.25\tsmall_variant\trs2\n')
if __name__=='__main__':unittest.main()
