import os
from pathlib import Path
from argparse import ArgumentParser

parser = ArgumentParser(description="Download Dataset")

parser.add_argument(
    "-i",
    "--url",
    default="",
    help="Indicate the url the dataset, e.g., http://konect.cc/networks/douban/"
)

def _main():
    args = parser.parse_args()
    graph_name = ""
    if args.url != "" and args.url[-1] == "/":
        graph_name = args.url[:-1].split('/')[-1]
    elif args.url != "" and args.url[-1] != "/":
        graph_name = args.url.split('/')[-1]
    else:
        raise ValueError("Please input a valid url")
    zip_url = "http://konect.cc/files/download.tsv.%s.tar.bz2"%graph_name
    zip_name = graph_name + ".tar.bz2"
    if os.path.exists(f"./{graph_name}/"):
        return
    if os.system(f"wget {zip_url} -O {zip_name} > /dev/null 2>&1") != 0:
        print(f"Fail to download the dataset {graph_name}!")
        return 
    if os.system(f"mkdir {graph_name} > /dev/null 2>&1") != 0:
        print(f"Fail to mkdir {graph_name}!")
        return 
    if os.system(f"tar -xjf {zip_name} > /dev/null 2>&1") != 0:
        print(f"Fail to decompress zip {zip_name}!")
        return 
    os.remove(zip_name)
    if not os.path.exists(f"{graph_name}/out.{graph_name}"):
        os.removedirs(graph_name)
        return     
    os.rename(f"{graph_name}/out.{graph_name}", f"{graph_name}/graph.txt")
    
    print(zip_url)

if __name__ == "__main__":
    _main()