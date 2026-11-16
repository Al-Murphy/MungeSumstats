"""OS process memory for benchmark safety, including compressed/swapped pages.

macOS layout follows sys/resource.h rusage_info_v0 and libproc.h. The physical
footprint includes compressed/private accounting that RSS alone omits. External
time remains authoritative for separately reported peak process RSS.
"""
import ctypes,os,pathlib,sys

class _UsageV0(ctypes.Structure):
    _fields_=[('uuid',ctypes.c_uint8*16)]+[(name,ctypes.c_uint64) for name in
        ('user_time','system_time','pkg_idle_wkups','interrupt_wkups','pageins','wired_size',
         'resident_size','phys_footprint','proc_start_abstime','proc_exit_abstime')]

_LIBPROC=None
if sys.platform=='darwin':
    _LIBPROC=ctypes.CDLL('/usr/lib/libproc.dylib',use_errno=True)
    _LIBPROC.proc_pid_rusage.argtypes=[ctypes.c_int,ctypes.c_int,ctypes.c_void_p]
    _LIBPROC.proc_pid_rusage.restype=ctypes.c_int

def process_memory(pid):
    if sys.platform=='darwin':
        usage=_UsageV0()
        if _LIBPROC.proc_pid_rusage(pid,0,ctypes.byref(usage)):
            raise OSError(ctypes.get_errno(),'proc_pid_rusage failed',pid)
        resident=int(usage.resident_size);footprint=int(usage.phys_footprint)
        return dict(resident_bytes=resident,physical_footprint_bytes=footprint,
                    guard_bytes=max(resident,footprint),scope='macOS max(RSS,phys_footprint including compressed/private accounting)')
    if sys.platform.startswith('linux'):
        values={}
        for line in pathlib.Path(f'/proc/{pid}/status').read_text().splitlines():
            if line.startswith(('VmRSS:','VmSwap:')):
                key,value,*_=line.split();values[key.rstrip(':')]=int(value)*1024
        resident=values.get('VmRSS',0);swapped=values.get('VmSwap',0)
        return dict(resident_bytes=resident,physical_footprint_bytes=None,
                    guard_bytes=resident+swapped,scope='Linux RSS+VmSwap')
    raise RuntimeError('Benchmark safety monitor supports macOS/Linux; runtime MSS does not require it')
