from __future__ import annotations
import errno,hashlib,json,os,shutil,signal,subprocess,sys,time
from pathlib import Path
from process_memory import process_memory
ROOT=Path(__file__).resolve().parent

def sha256(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()

def pin(path):
 p=Path(path).resolve();st=p.stat();h=hashlib.sha256();md=hashlib.md5();lines=0
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b);md.update(b);lines+=b.count(b'\n')
 if (p.stat().st_size,p.stat().st_mtime_ns)!=(st.st_size,st.st_mtime_ns):raise RuntimeError('Pinned file changed while hashing')
 return dict(path=str(p),bytes=st.st_size,mtime_ns=st.st_mtime_ns,sha256=h.hexdigest(),md5=md.hexdigest(),newline_count=lines)

def write_once(path,value):
 with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')

def check_pin(info):
 p=Path(info['path'])
 if p.stat().st_size!=info['bytes'] or sha256(p)!=info['sha256']:raise RuntimeError('Immutable file pin changed: '+str(p))

def disk_bytes(root):
 total=0
 for directory,dirs,files in os.walk(root,followlinks=False):
  for name in files:
   try:total+=(Path(directory)/name).lstat().st_size
   except FileNotFoundError:pass
 return total

def tree_usage(pid):
 lines=subprocess.check_output(['ps','-axo','pid=,ppid='],text=True)
 rows=[tuple(map(int,line.split())) for line in lines.splitlines() if line.strip()]
 live={p for p,_ in rows};descendants={pid}
 while True:
  more=descendants|{p for p,parent in rows if parent in descendants}
  if more==descendants:break
  descendants=more
 rss=footprint=guard=0
 for child in descendants:
  try:use=process_memory(child)
  except OSError as exc:
   if exc.errno==errno.ESRCH or child not in live:continue
   raise
  rss+=use['resident_bytes'];footprint+=use['physical_footprint_bytes'] or 0;guard+=use['guard_bytes']
 return dict(rss_bytes=rss,physical_footprint_bytes=footprint,conservative_guard_bytes=guard,processes=len(descendants))

def guarded(command,case,receipt,log,config,deadline,env=None):
 if Path(receipt).exists() or Path(log).exists():raise FileExistsError('Existing guard/log refused')
 if shutil.disk_usage(case).free<config['min_filesystem_free_bytes']:raise RuntimeError('160GiB reserve unavailable')
 # Fail before study processing if this operator cannot measure the full tree.
 tree_usage(os.getpid())
 start=time.monotonic();reason=None;peaks=dict(rss_bytes=0,physical_footprint_bytes=0,conservative_guard_bytes=0,processes=0);samples=0;minimum_free=shutil.disk_usage(case).free;peak_disk=0
 with Path(log).open('x') as out:
  child=subprocess.Popen(command,stdout=out,stderr=subprocess.STDOUT,start_new_session=True,env=env)
  try:
   while child.poll() is None:
    use=tree_usage(child.pid);samples+=1
    for key in peaks:peaks[key]=max(peaks[key],use[key])
    free=shutil.disk_usage(case).free;minimum_free=min(minimum_free,free);growth=disk_bytes(case);peak_disk=max(peak_disk,growth)
    if any(use[k]>config['max_tree_memory_bytes'] for k in ['rss_bytes','physical_footprint_bytes','conservative_guard_bytes']):reason='8GiB process-tree memory target exceeded'
    elif free<config['min_filesystem_free_bytes']:reason='160GiB filesystem free reserve reached'
    elif growth>config['max_case_disk_bytes']:reason='Owned case disk-growth cap reached'
    elif time.monotonic()>deadline:reason='Owned case deadline reached'
    if reason:
     if child.poll() is None:os.killpg(child.pid,signal.SIGKILL)
     break
    time.sleep(config['poll_seconds'])
   code=child.wait()
  except BaseException as exc:
   if child.poll() is None:os.killpg(child.pid,signal.SIGKILL);child.wait()
   code=child.returncode;reason='Guard measurement/operator failure: '+type(exc).__name__+': '+str(exc)
 state=dict(command=command,pid=child.pid,status='censored' if reason else ('success' if code==0 else 'failed'),exit_code=code,reason=reason,natural_oom=False,wall_seconds=time.monotonic()-start,samples=samples,peak_process_tree=peaks,min_free_bytes=minimum_free,peak_owned_case_bytes=peak_disk,limits=dict(memory_bytes=config['max_tree_memory_bytes'],minimum_free_bytes=config['min_filesystem_free_bytes'],case_disk_bytes=config['max_case_disk_bytes']),scope='Full owned process tree: RSS, macOS physical footprint and sum(max(RSS,footprint)); Linux guard includes swap. No other workload is signaled.')
 write_once(receipt,state)
 if reason or code:raise RuntimeError('Guarded command did not complete: '+str(state))
 return state
