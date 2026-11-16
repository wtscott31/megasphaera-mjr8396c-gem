from common import *
import concurrent.futures,sys
def execute(pair):return run([sys.executable,str(ROOT/'scripts/matched_micom_v6.py'),'--pair',pair],'matched_micom_'+pair,timeout=7200)
if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows=list(pool.map(execute,[n+'_matched_curated_sensitivity' for n in ['MJR8396C','NM10','BL7']]))
    record('results/matched_micom_v6/curated_parallel_execution.json',rows)
