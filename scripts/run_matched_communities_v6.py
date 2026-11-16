from common import *
import concurrent.futures,sys
def execute(pair):
    cmd=[sys.executable,str(ROOT/'scripts/matched_micom_v6.py'),'--pair',pair]
    return run(cmd,'matched_micom_'+pair,timeout=7200)
if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        rows=list(pool.map(execute,['MJR8396C','NM10','BL7','GPRcurated_original22L']))
    record('results/matched_micom_v6/parallel_execution.json',rows)
