import importlib.util,json,os,pathlib,signal,subprocess,sys,time,unittest

path=pathlib.Path(__file__).resolve().parents[2]/'scripts/process_memory.py'
spec=importlib.util.spec_from_file_location('process_memory',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class ActualAllocationGuard(unittest.TestCase):
    def test_small_actual_child_is_preserved_then_killed_above_limit(self):
        child_code='''import sys,time
x=bytearray(16*1024**2)
for i in range(0,len(x),4096):x[i]=1
print('READY',flush=True)
sys.stdin.readline()
y=bytearray(96*1024**2)
for i in range(0,len(y),4096):y[i]=1
time.sleep(30)
'''
        child=subprocess.Popen([sys.executable,'-c',child_code],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
        try:
            self.assertEqual(child.stdout.readline().strip(),'READY')
            limit=64*1024**2
            small=module.process_memory(child.pid)
            self.assertLess(small['guard_bytes'],limit)
            self.assertIsNone(child.poll())
            child.stdin.write('allocate\n');child.stdin.flush()
            observed=None
            for _ in range(200):
                observed=module.process_memory(child.pid)
                if observed['guard_bytes']>limit:
                    os.kill(child.pid,signal.SIGKILL);break
                time.sleep(.02)
            self.assertIsNotNone(observed)
            self.assertGreater(observed['guard_bytes'],limit)
            self.assertEqual(child.wait(timeout=5),-signal.SIGKILL)
            print(json.dumps(dict(control=small,trigger=observed,limit_bytes=limit,exitcode=child.returncode)))
        finally:
            if child.poll() is None:child.kill();child.wait()
            child.stdin.close();child.stdout.close()

if __name__=='__main__':unittest.main()
