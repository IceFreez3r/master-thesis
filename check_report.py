import os
import sys
import glob


results_dir = "/project/hfa_work/ENCODE/code/snakemake-pipeline/results/sqanti"
TOOLS = ["flair", "isoquant", "isotools", "isotools_new_tss", "stringtie"]
TISSUES = ["aorta", "brain", "colon", "heart", "lung", "muscle"]

def get_tools(argv):
    if len(argv) >= 1 and argv[0] in TOOLS:
        return [argv[0]]
    elif len(argv) >= 2 and argv[1] in TOOLS:
        return [argv[1]]
    else:
        return TOOLS

def get_tissues(argv):
    if len(argv) >= 1 and argv[0] in TISSUES:
        return [argv[0]]
    elif len(argv) >= 2 and argv[1] in TISSUES:
        return [argv[1]]
    else:
        return TISSUES

if len(sys.argv) < 2 or sys.argv[1] == "list":
    for root, dirs, files in os.walk(results_dir):
        for file in files:
            if file.endswith(".html"):
                print(os.path.join(root, file))
elif sys.argv[1] == "check":
    tools = get_tools(sys.argv[2:])
    tissues = get_tissues(sys.argv[2:])
    for tissue in tissues:
        for tool in tools:
            report_path = os.path.join(results_dir, tool, "qc", tissue, "*.html")
            for file in glob.glob(report_path):
                print(tissue, tool, file, sep='\t')
elif sys.argv[1] == "open":
    tools = get_tools(sys.argv[2:])
    tissues = get_tissues(sys.argv[2:])
    cmd = "firefox"
    for tissue in tissues:
        for tool in tools:
            report_path = os.path.join(results_dir, tool, "qc", tissue, "*.html")
            for file in glob.glob(report_path):
                cmd += f" --new-tab {file}"
    os.system(cmd)
elif sys.argv[1] == "copy":
    if len(sys.argv) < 3:
        print("Usage: python check_report.py copy <destination> [tool] [tissue]")
        sys.exit
    destination = sys.argv[2]
    os.makedirs(destination, exist_ok=True)

    tools = get_tools(sys.argv[3:])
    tissues = get_tissues(sys.argv[3:])
    for tissue in tissues:
        for tool in tools:
            classification_path = os.path.join(results_dir, tool, "qc", tissue, "*_classification.txt")
            for file in glob.glob(classification_path):
                os.system(f"cp {file} {os.path.join(destination, tissue + '_' + tool + '_' + os.path.basename(file).split('_')[1])}")
