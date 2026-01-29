# This is a script to download data using SRA accession numbers, see the main() function
# at the bottom for what arguments you need to edit

#import necessary libraries
import re
import subprocess
import os
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
import threading

class SraDownloader:
    def __init__(self, out_dir: str, temp_dir: str, sra_dir: str, threads: int = 4, dry_run: bool = False):
        self.dry_run = dry_run
        self.out_dir = out_dir
        self.temp_dir = temp_dir
        self.sra_dir = sra_dir
        self.threads = threads

    def prefetch_one(self, acc: str):
        cmd = ["prefetch", acc, "--output-directory", self.sra_dir]
        run_command(cmd, dry_run=self.dry_run)
        print(f"prefetch successful for {acc}")

    def fasterq_dump_one(self, acc: str):
        out_dir_acc = f"{self.out_dir}/{acc}"
        tmp_dir_acc = f"{self.temp_dir}/{acc}"
        sra_path = f"{self.sra_dir}/{acc}/{acc}.sra"

        os.makedirs(out_dir_acc, exist_ok=True)
        os.makedirs(tmp_dir_acc, exist_ok=True)

        # In dry-run, the .sra won't exist yet; only enforce this on real runs.
        if not self.dry_run and not os.path.exists(sra_path):
            raise FileNotFoundError(f"Missing SRA file: {sra_path}")

        cmd = [
            "fasterq-dump", sra_path,
            "--outdir", out_dir_acc,
            "--temp", tmp_dir_acc,
            "--threads", str(self.threads),
            "--split-files"
        ]
        run_command(cmd, dry_run=self.dry_run)
        print(f"fasterq-dump successful for {acc}")

    def is_done(self, acc: str) -> bool:
        out_dir_acc = f"{self.out_dir}/{acc}"
        se = f"{out_dir_acc}/{acc}.fastq"
        pe1 = f"{out_dir_acc}/{acc}_1.fastq"
        if os.path.exists(se) and os.path.getsize(se) > 0:
            return True
        elif os.path.exists(pe1) and os.path.getsize(pe1) > 0:
            return True
        else:
            return False

#logs
LOG_LOCK = threading.Lock()

def append_log(log_path: str, acc: str, status: str, msg: str):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"{ts}\t{acc}\t{status}\t{msg}\n"
    # Ensure parent directory exists
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    # append logs
    with LOG_LOCK:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(line)

def read_accessions(path: str) -> list[str]:
    acc_list = []
    with open(path, "r") as file:
        for line in file:
            #handle inline comments
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            # Split on any whitespace
            tokens = line.split()
            acc_list.extend(tokens)

    return acc_list

def validate_accessions(acc: str) -> bool:
    return re.fullmatch(r"(SRR|ERR|DRR)\d+", acc) is not None


def dedup_preserve_order(items: list[str]) -> list[str]:
    seen = set()
    out = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out

def process_one_accession(acc: str, out_dir: str, tmp_dir: str, sra_dir: str, threads: int,
                          dry_run: bool, delete_tmp_after: bool, delete_sra_after: bool, log_path: str) -> tuple[str, bool, str]:
    downloader = SraDownloader(out_dir, tmp_dir, sra_dir, threads=threads, dry_run=dry_run)

    # Skip if already done
    if downloader.is_done(acc):
        append_log(log_path, acc, "SKIP", "fastq already downloaded")
        return acc, True, "SKIP (fastq already downloaded)"

    try:
        downloader.prefetch_one(acc)
        downloader.fasterq_dump_one(acc)

        tmp_dir_acc = f"{tmp_dir}/{acc}"
        sra_dir_acc = f"{sra_dir}/{acc}"

        if delete_tmp_after:
            safe_rmtree(tmp_dir_acc, dry_run=dry_run)

        if delete_sra_after:
            safe_rmtree(sra_dir_acc, dry_run=dry_run)

        append_log(log_path, acc, "OK", "prefetch + fasterq-dump completed")
        return acc, True, "OK"

    except Exception as e:
        append_log(log_path, acc, "FAIL", str(e).replace("\n", " "))
        return acc, False, str(e)

def safe_rmtree(path: str, dry_run: bool = False):
    if not os.path.exists(path):
        return
    print(f"Cleaning: {path}")
    if dry_run:
        return
    shutil.rmtree(path)

def run_command(cmd: list[str], dry_run: bool = False):
    printable = " ".join(cmd)
    print(f"Running: {printable}")
    if dry_run:
        return None

    result = subprocess.run(
        cmd,
        text=True,
        capture_output=True
    )

    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {printable}\n{result.stderr}")

    return result

def main():
    #manually choose the paths you want (accessions, output, temp & sra storage, logs)
    test_path = "/mnt/c/Users/eddie/Personal_Research/fastq_accessions/ash1_RNA.txt"
    out_dir = "/mnt/d/Personal_Research/Data/fastq_files"
    tmp_dir = "/mnt/d/Personal_Research/Data/tmp_files"
    sra_dir = "/mnt/d/Personal_Research/Data/sra_cache"
    log_path = "/mnt/c/Users/eddie/Personal_Research/logs/sra_fetch.log"

    #choose whether you want sra/tmp files deleted after job completion
    delete_tmp_after = True
    delete_sra_after = True

    # choose concurrency + threads here. If you are running on a local machine, ensure that
    # it can handle the number of concurrent jobs (jobs * threads <= cpu cores).
    # Also choose if you want a dry test run
    jobs = 2
    threads = 2
    dry_run = False

    #read in desired accessions here
    acc_list = read_accessions(test_path)
    new_list = []
    ok = 0
    failed = 0

    #remove invalid accessions
    for acc in acc_list:
        acc = acc.strip().upper()
        if not validate_accessions(acc):
            print(f"{acc} is not a valid accession, removed from list")
            continue
        else:
            print(f"{acc} validated")
            new_list.append(acc)

    #Remove duplicate accessions
    acc_list = dedup_preserve_order(new_list)

    #if accession list is now empty, stop
    if not acc_list:
        print("No accessions found")
        return

    #Run concurrent downloads (don't change anything here)
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futures = [
            ex.submit(process_one_accession, acc, out_dir, tmp_dir, sra_dir, threads, dry_run, delete_tmp_after, delete_sra_after, log_path)
            for acc in acc_list
        ]
        for fut in as_completed(futures):
            acc, success, msg = fut.result()
            if success:
                ok += 1
                print(f"{acc}: {msg}")
            else:
                failed += 1
                print(f"{acc}: FAILED - {msg}")
    print(f"Done. Successful: {ok}, Failed: {failed}, Total: {len(acc_list)}")

if __name__ == "__main__": main()
