import os
import re
import json
import tempfile
import logging
import matplotlib.pyplot as plt
from argparse import ArgumentParser

OUTPUT_PATH = "./out/"
TOOL_PATH = "./src/btfc"
# DATASET = ["bi-lj", "bag-nytimes", "deli-ui", "lastfm_band", "orkut", "yahoo"]
# DATASET = ["orkut-groupmemberships"]
DATASET = ["yahoo-song"]
NUM_OF_DUMPLICATE = 20

parser = ArgumentParser(description="Butterfly Counter Running Script")

parser.add_argument(
    "-m",
    "--method",
    default="OUR",
    help="Indicate the name of method."
)

parser.add_argument(
    "-e",
    "--experiment",
    default="s1s2",
    help="Indicate the name of experiments."
)

class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"

def warn(message):
    print(bcolors.WARNING + "[WARN] " + message + bcolors.ENDC)

def error(message):
    print(bcolors.FAIL + "[ERROR] " + message + bcolors.ENDC)

def info(message):
    print(bcolors.HEADER + "[INFO] " + message + bcolors.ENDC)

def ok(message):
    print(bcolors.OKGREEN + "[OK] " + message + bcolors.ENDC)

def fail(message):
    print(bcolors.FAIL + "[FAIL] " + message + bcolors.ENDC)

def median(lst):
    if lst == []:
        return 0

    n = len(lst)
    s = sorted(lst)
    return (s[n//2-1]/2.0+s[n//2]/2.0, s[n//2])[n % 2] if n else None

# def get_ground_truth(benchmark, d):
#     if d == 10:
#         for filename, gt in GT_dict.items():
#             if benchmark in filename:
#                 return gt
#     else:
#         cmd = "./src/btfc -m exact ./dataset/{} -d {} > gt.txt".format(benchmark, d)
#         # print(cmd)
#         os.system(cmd)
#         with open("gt.txt", 'r') as file:
#             for line in file:
#                 match = re.search(r'BTF: (\d+)', line)
#                 if match:
#                     number = int(match.group(1))
#                     return number
#         return 0

def generate_ground_truth(density=10):
    GROUND_TRUTH_PATH = OUTPUT_PATH + "ground-truth.json"
    ground_truth = dict()

    if os.path.exists(GROUND_TRUTH_PATH):
        with open(OUTPUT_PATH + "ground-truth.json", "r") as f:
            try:
                ground_truth = json.load(f)
            except:
                warn("Unsupported json format in ground-truth.json")
                ground_truth = dict()
    ## find the exact butterfly number for each benchmark
    for benchmark in DATASET:
        info("Generating ground truth for {}..".format(benchmark))
        if benchmark not in ground_truth.keys() or ground_truth[benchmark] == 0:
            if not os.path.exists(f"./dataset/{benchmark}"):
                error(f"Please first download the graph file of {benchmark}")
            if not os.path.exists(f"./dataset/{benchmark}/graph-sort.bin"):
                warn(f"No graph-sort.bin found in ./dataset/{benchmark}/")
                cmd = "{} -m txt2bin -p ./dataset/{}/ -d {} > {}".format(
                TOOL_PATH, benchmark, density, output_path)
                os.system(cmd)
                ok(f"Convert txt to bin for {benchmark} successfully!")
            output_path = "{}{}_gt.log".format(OUTPUT_PATH, benchmark)
            cmd = "{} -m exact -p ./dataset/{}/ -d {} > {}".format(
                TOOL_PATH, benchmark, density, output_path)
            print(cmd)
            os.system(cmd)
            with open(output_path, 'r') as file:
                for line in file:
                    match = re.search(r'BTF: (\d+)', line)
                    if match:
                        number = int(match.group(1))
                        ground_truth[benchmark] = number
                        ok("Ground truth of {} is {}".format(benchmark, number))
                        break
    with open(OUTPUT_PATH + "ground-truth.json", "w") as f:
        json.dump(ground_truth, f)

# run the two-leveled sampling method
def run_tls(benchmark, d, s1, s2, rbase):
    if not os.path.exists("./out/tls"):
        os.mkdir("./out/tls")
    if not os.path.exists("./out/tls/{}".format(benchmark)):
        os.mkdir("./out/tls/{}".format(benchmark))
    if not os.path.exists(f"./dataset/{benchmark}/graph-sort.bin"):
        fail(f"Fail to find benchmark {benchmark}")
        raise ValueError(f"Fail to find graph-sort-bin for {benchmark}")
    log_file_path = "./out/tls/{}/{}sqrtm_{}sqrtm.log".format(
        benchmark, "" if s1 == 1 else s1, "" if s2 == 1 else s2)
    nums = ' '.join(str(i) for i in range(1, NUM_OF_DUMPLICATE))
    os.system("for i in {}; do ./src/btfc -m tls -p ./dataset/{} -d {} --s1 {} --s2 {} --rbase {} -q; done > {}".format(
        nums, benchmark, d, s1, s2, rbase, log_file_path))
    with open(log_file_path, "r") as f:
        log_str = "".join(f.readlines())
        BTF = re.findall("BTF: (.+?),", log_str)
        Time = re.findall("Time: (.+?\..+?),", log_str)
        Query = re.findall("Query: (.+?),", log_str)
        
        BTF, Time, Query = list(map(float, BTF)), list(map(float, Time)), list(map(float, Query))
    
    return {"BTF": BTF, "Time": Time, "Query": Query}

def run_wps(benchmark, d, rround=10):
    if not os.path.exists("./out/wps"):
        os.mkdir("./out/wps")
    if not os.path.exists("./out/wps/{}".format(benchmark)):
        os.mkdir("./out/wps/{}".format(benchmark))
    if not os.path.exists(f"./dataset/{benchmark}/graph-sort.bin"):
        fail(f"Fail to find benchmark {benchmark}")
        raise ValueError(f"Fail to find graph-sort-bin for {benchmark}")
    log_file_path = "./out/wps/{}/{}.log".format(benchmark, benchmark)
    nums = ' '.join(str(i) for i in range(1, NUM_OF_DUMPLICATE))
    os.system("for i in {}; do ./src/btfc -m wps -p ./dataset/{} -d {} --rround {} -q ; done > {}".format(
        nums, benchmark, d, rround, log_file_path))
    with open(log_file_path, "r") as f:
        log_str = "".join(f.readlines())
        BTF = re.findall("BTF: (.+?),", log_str)
        Time = re.findall("Time: (.+?\..+?),", log_str)
        Query = re.findall("Query: (.+?),", log_str)
        
        BTF, Time, Query = list(map(float, BTF)), list(map(float, Time)), list(map(float, Query))
    
    return {"BTF": BTF, "Time": Time, "Query": Query}

def output_statics(benchmark, method, performance):
    GROUND_TRUTH_PATH = OUTPUT_PATH + "ground-truth.json"
    ground_truth = dict()
    with open(GROUND_TRUTH_PATH, "r") as f:
        ground_truth = json.load(f)
    
    if benchmark not in ground_truth.keys():
        error(f"Unknown ground truth for {benchmark}")

    btf_error_rates = [(btf - ground_truth[benchmark]) / ground_truth[benchmark] for btf in performance["BTF"]]

    median_btf = median(performance["BTF"])
    btf_error_rate = abs(median_btf - ground_truth[benchmark])/ground_truth[benchmark]
    median_time = median(performance["Time"])
    median_query = median(performance["Query"])
    
    ok("{}: {}, {}, {}, {}".format(benchmark, method, btf_error_rate, median_time, median_query))

    return btf_error_rates
    
def _main():
    args = parser.parse_args()

    if not os.path.exists(OUTPUT_PATH):
        os.mkdir(OUTPUT_PATH)
    
    if args.experiment == "ground-truth":
        generate_ground_truth()
    elif args.experiment == "comp-all":
        for benchmark in DATASET:
            # try:
            tls_performance = run_tls(benchmark, 10, 0.2, 1, 100)
            y1 = output_statics(benchmark, "TLS", tls_performance)
            wps_performance = run_wps(benchmark, 10, rround=10)
            y2 = output_statics(benchmark, "WPS", wps_performance)

            # tls_mean_error = sum(tls_performance["BTF"])/len(tls_performance["BTF"])
            # tls_variance = sum((x - tls_mean_error) ** 2 for x in tls_performance["BTF"])
            # wps_mean_error = sum(wps_performance["BTF"])/len(wps_performance["BTF"])
            # wps_variance = sum((x - wps_mean_error) ** 2 for x in wps_performance["BTF"])

            # print(f"TLS: {tls_mean_error}, {tls_variance}")
            # print(f"WPS: {wps_mean_error}, {wps_variance}")

            # x = list(range(1, NUM_OF_DUMPLICATE))
            x = [1 for _ in range(NUM_OF_DUMPLICATE -1 )]
            # y1 = btf_error_rates
            # y1 = tls_performance["BTF"]
            # y2 = wps_performance["BTF"]
            plt.scatter(x, y1, color='blue', label='TLS', marker='o')
            plt.scatter(x, y2, color='red', label='WPS', marker='x')

            plt.xlabel('# of run')
            plt.ylabel('Relative Error')
            plt.title('Comparison of TLS and WPS')
            plt.legend()
            plt.grid(True)
            plt.show()
            # plt.savefig(f"./out/{benchmark}_comp.png")

            # except:
            #     fail(f"Fail to run estimate methods on {benchmark}")
            
if __name__ == "__main__":
    _main()