"""
Ensure consistent directory names

"""

import attrmap as ap
import os

dir = ap.AttrMap()

# Workflow dirs
dir.env = os.path.join(workflow.basedir, "envs")
dir.rules = os.path.join(workflow.basedir, "rules")
dir.scripts = os.path.join(workflow.basedir, "scripts")

# Output locations
dir.out.base = config.args.output
dir.out.json = os.path.join(dir.out.base, "json")
dir.out.download = os.path.join(dir.out.base, "download")

# Logs
dir.logs = os.path.join(dir.out.base, "logs")
# versions
dir.versions = os.path.join(dir.out.base, "versions")
