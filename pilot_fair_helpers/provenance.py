import json
import subprocess
import sys
from pathlib import Path
from datetime import datetime

def run(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode().strip()
    except subprocess.CalledProcessError as e:
        return f"ERROR: {e.output.decode().strip()}"

def main():
    wf = sys.argv[1]
    wf_path = Path(wf).resolve()

    repo = run(f"git -C '{wf_path.parent}' rev-parse --show-toplevel")

    git_info = {
        "timestamp": datetime.utcnow().isoformat(),
        "repo": repo,
        "workflow": wf_path.name,
        "remote": run(f"git -C '{repo}' remote get-url origin"),
        "commit": run(f"git -C '{repo}' rev-parse HEAD"),
        "version": run(f"git -C '{repo}' describe --tags --always")
    }

    # write git info to JSON file
    with open("git_info.json", "w") as f:
        json.dump(git_info, f, indent=2)

    # generate CWL snapshot
    subprocess.run(
        f"cwltool --pack '{wf}' > workflow_snapshot.json",
        shell=True
    )

    print("Wrote git_info.json and workflow_snapshot.json")

if __name__ == "__main__":
    main()