import subprocess
import re
import time
import psutil
import os

def get_resource_output():

    # Run the `top` command and capture its output
    command = '''free | awk '/Mem:/ {print $3/$2 * 100.0}' && top -bn1 | grep "Cpu(s)" | sed "s/.*, *\\([0-9.]*\\)%* id.*/\\1/" | awk '{print 100 - $1}' '''
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, text=True)
    mem_usage, cpu_usage = result.stdout.splitlines()
    return (float(mem_usage), float(cpu_usage))

# def parse_top_output(output):
#     lines = output.splitlines()
#     matrix = []
#     header = None
#     # Find the header line and the data lines
#     for line in lines:
#         if line.startswith("PID"):
#             header = re.split(r'\s+', line.strip())
#         elif re.match(r'^\s*\d+', line):
#             data = re.split(r'\s+', line.strip())
#             matrix.append(data)
#     return header, matrix

# def get_top_processes(matrix, top_n):
#     # Sort processes by %CPU usage (index 8) in descending order and take the top N
#     sorted_processes = sorted(matrix, key=lambda x: float(x[8]), reverse=True)[:top_n]
#     return sorted_processes

# def compute_totals(processes):
#     cpu_total = sum(float(p[8]) for p in processes) / psutil.cpu_count()
#     mem_total = sum(float(p[9]) for p in processes)
#     return cpu_total, mem_total

def log_totals(cpu_total, mem_total, path, timestamp, ttime):
    with open(path+"cpu_memory_"+str(timestamp)+".txt", 'a') as f:  # Append mode
        f.write(f"TIME:{ttime}|CPU:{cpu_total:.4f}|MEM:{mem_total:.4f}\n")

def TOP_LOGGER(interval, path, timestamp):
    
    # Set some finite duration so it doesn't infinitely leak memory
    end_time = time.time() + 100000

    # Initialize time
    ttime = 0

    # Check if output file already exists, if so remove
    if os.path.exists(path+"cpu_memory_"+str(timestamp)+".txt"):
        os.remove(path+"cpu_memory_"+str(timestamp)+".txt")

    # Run logger
    while time.time() < end_time:
        mem_usage, cpu_usage = get_resource_output()
        log_totals(cpu_usage, mem_usage, path, timestamp, ttime)
        time.sleep(interval)
        ttime += interval

if __name__ == "__main__":
    TOP_LOGGER(2, "output/", 10000)