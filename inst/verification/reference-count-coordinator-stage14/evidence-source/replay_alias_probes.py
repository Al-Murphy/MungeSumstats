"""Run byte-preserved review probes with only ROOT/CLONE literals relocated."""
import ast
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
OLD_ROOT='/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-fixes-independent-stage11'
OLD_CLONE='/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-fixes-stage11'


def main():
    name,phase=sys.argv[1:]
    assert name in ('probe_alias_gate.py','probe_guard_order_and_data_link.py')
    assert phase in ('red','green')
    source=HERE/'alias-evidence/review-probes'/name
    target=HERE/'alias-evidence'/phase/('alias' if name=='probe_alias_gate.py' else 'order')
    target.mkdir(parents=True)
    mapping={OLD_ROOT:str(target),OLD_CLONE:str(HERE)}
    class Relocate(ast.NodeTransformer):
        def visit_Constant(self,node):
            if isinstance(node.value,str) and node.value in mapping:
                return ast.copy_location(ast.Constant(mapping[node.value]),node)
            return node
    tree=Relocate().visit(ast.parse(source.read_bytes(),filename=str(source)))
    ast.fix_missing_locations(tree)
    record={'probe_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'path_only_literal_replacements':mapping,'probe_logic_changed':False,'phase':phase}
    (target/'relocation.json').write_text(json.dumps(record,indent=2)+'\n')
    exec(compile(tree,str(source),'exec'),{'__name__':'__main__','__file__':str(source)})


if __name__=='__main__':main()
