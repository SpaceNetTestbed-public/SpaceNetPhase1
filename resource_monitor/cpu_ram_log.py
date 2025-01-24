import psutil
import argparse
import time
import csv
import os
import datetime

def bytesto(bytes, to, bsize=1024):
    """Convert bytes to megabytes, etc."""
    a = {'k': 1, 'm': 2, 'g': 3, 't': 4, 'p': 5, 'e': 6}
    r = float(bytes)
    for _ in range(a[to]):
        r /= bsize
    return r

def get_process_info(proc, username):
    """Retrieve and filter process information."""
    try:
        process_info = proc.as_dict(attrs=["pid", "name", "username", "memory_info"])
        if process_info["name"] == "sleep":
            return None
        if process_info["username"] == username and process_info["name"] == "python3" and proc.memory_info()[0] > 50e+2:
            memory_percent = proc.memory_info()[0]
            cpu_percent = proc.cpu_percent(1)
            return memory_percent, cpu_percent
    except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
        pass
    return None

def log_resource_usage(args):
    """Log CPU and memory usage to a file."""
    current_count = 0
    header = ['datetime', 'CPU [%]', 'Memory [%]']

    try:
        print("Start data logging")

        while current_count < args.logcount:
            memory_avg = []
            cpu_avg = []

            for _ in range(args.interval):
                memory = []
                cpu_usage = []

                for proc in psutil.process_iter():
                    result = get_process_info(proc, args.username)
                    if result:
                        memory_percent, cpu_percent = result
                        memory.append(memory_percent)
                        cpu_usage.append(cpu_percent)

                memory_avg.append(sum(memory))
                cpu_avg.append(sum(cpu_usage))
                time.sleep(1)

            final_cpu = round(sum(cpu_avg) / len(cpu_avg), 4)
            final_memory = round(sum(memory_avg) / len(memory_avg), 4) / psutil.virtual_memory().total * 100

            with open(args.output, "a", newline="", encoding="UTF-8") as csv_file:
                writer = csv.writer(csv_file, delimiter='\t')
                if os.path.getsize(args.output) == 0:
                    writer.writerow(header)
                writer.writerow([datetime.datetime.now().strftime("%H:%M:%S %d/%m/%Y"), final_cpu, final_memory])

            current_count += 1

    except KeyboardInterrupt:
        print(f"Logging finished and saved to file {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--username", type=str, required=True,
                        help="Specify the username for which we record CPU and Memory usage")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to the file where the CPU and Memory usage gets logged (in tsv format)")
    parser.add_argument("--interval", type=int, required=True,
                        help="Specify the time (in seconds) over which the CPU and Memory usage gets averaged")
    parser.add_argument("--logcount", type=int, required=True,
                        help="Specify the number of logs")
    args = parser.parse_args()

    log_resource_usage(args)