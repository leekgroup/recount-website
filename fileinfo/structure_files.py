#!/usr/bin/env python
"""
structure_files.py

Symlinks to create directory structure for syncing Recount2 files with Amazon
Cloud Drive. Requires upload_table.tsv (in current directory).

Arg 1: root dir for files to upload
Arg 2: location of output of junctions_by_project.py and add_knowngene.py
"""
import sys
import os
from glob import glob

root_dir = sys.argv[1]
junctions_path = sys.argv[2]

# Create root dir if it doesn't exist
try:
    os.makedirs(root_dir)
except OSError:
    if os.path.exists(root_dir):
        pass
    else:
        raise

current_project = None
for line in sys.stdin:
    filename, basename, project = line.strip().split('\t')
    if project != current_project:
        project_dir = os.path.join(root_dir, project)
        bw_dir = os.path.join(root_dir, project, 'bw')
        os.makedirs(bw_dir)
        current_project = project
    if basename.endswith('.bw'):
        os.link(filename, os.path.join(bw_dir, basename))
    else:
        os.link(filename, os.path.join(project_dir, basename))

os.link(
        os.path.join(junctions_path, 'sample_ids.tsv'),
        os.path.join(root_dir, 'sample_ids.tsv')
    )
