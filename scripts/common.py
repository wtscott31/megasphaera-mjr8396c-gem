from pathlib import Path
import datetime, hashlib, json, subprocess, time, urllib.request

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'inputs/submission_reference'
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def record(name, value):
    p=ROOT/name; p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,indent=2,default=str)); return value
def fetch(url, rel, timeout=60):
    p=ROOT/rel; p.parent.mkdir(parents=True,exist_ok=True)
    row={'url':url,'path':rel,'retrieved_utc':now()}
    try:
        # System curl uses this Mac's trusted certificate store. Python's CA
        # bundle rejected the local proxy certificate; verification stays on.
        response=subprocess.run(['curl','--http1.1','-L','--fail','--silent','--show-error','--max-time',str(timeout),'-D',str(p)+'.headers',url,'-o',str(p)],capture_output=True,text=True,timeout=timeout+5)
        row['command']=['curl','--http1.1','-L','--fail','--max-time',str(timeout),url,'-o',rel]
        if response.returncode:raise RuntimeError(response.stderr.strip())
        row.update(bytes=p.stat().st_size,sha256=sha(p),outcome='downloaded')
    except Exception as e:
        row.update(outcome='failed',error=repr(e))
    with (ROOT/'logs/downloads.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(row['outcome'],url,flush=True)
    return p if row['outcome']=='downloaded' else None
def run(cmd,name,timeout=300,cwd=None):
    log=ROOT/'logs'/f'{name}.log';log.parent.mkdir(parents=True,exist_ok=True)
    row={'command':list(map(str,cmd)),'cwd':str(cwd or ROOT),'started_utc':now(),'log':str(log.relative_to(ROOT))}
    start=time.monotonic()
    try:
        with log.open('w') as f:
            f.write(json.dumps(row)+'\n');f.flush()
            r=subprocess.run(list(map(str,cmd)),cwd=cwd or ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
        row['returncode']=r.returncode
    except Exception as e:row['error']=repr(e)
    row.update(finished_utc=now(),seconds=time.monotonic()-start)
    with (ROOT/'logs/commands.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
    print(name,row.get('returncode',row.get('error')),flush=True)
    return row
