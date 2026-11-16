"""Replay exact original probe bytes with only the two documented paths relocated."""
import ast
import hashlib
import json
from pathlib import Path

ROOT=Path('/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-alias-fix-independent-stage11')
CLONE=Path('/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-alias-fix-stage11')
OLD_ROOT='/private/tmp/gwas-harvest-stage2/results/reference-count-coordinator-fixes-independent-stage11'
OLD_CLONE='/private/tmp/gwas-harvest-stage2/review-reference-count-coordinator-fixes-stage11'

for name,label in [('probe_alias_gate.py','alias'),('probe_guard_order_and_data_link.py','order')]:
    raw=(Path(OLD_ROOT)/name).read_bytes()
    copy=ROOT/name;copy.write_bytes(raw)
    target=ROOT/('independent-'+label);target.mkdir()
    mapping={OLD_ROOT:str(target),OLD_CLONE:str(CLONE)}
    class Relocate(ast.NodeTransformer):
        def __init__(self,mapping):self.mapping=mapping;self.count=0
        def visit_Constant(self,node):
            if isinstance(node.value,str) and node.value in self.mapping:
                self.count+=1
                return ast.copy_location(ast.Constant(self.mapping[node.value]),node)
            return node
    original=ast.parse(raw,filename=str(copy));transformer=Relocate(mapping)
    relocated=transformer.visit(ast.parse(raw,filename=str(copy)));ast.fix_missing_locations(relocated)
    assert transformer.count==2
    inverse=Relocate({v:k for k,v in mapping.items()}).visit(ast.parse(ast.unparse(relocated)))
    assert ast.dump(inverse,include_attributes=False)==ast.dump(original,include_attributes=False)
    record={'probe_sha256':hashlib.sha256(raw).hexdigest(),'literal_replacements':mapping,
            'changed_path_constants':2,'inverse_AST_roundtrip_exact':True,'probe_logic_changed':False}
    (target/'relocation.json').write_text(json.dumps(record,indent=2)+'\n')
    exec(compile(relocated,str(copy),'exec'),{'__name__':'__main__','__file__':str(copy)})
