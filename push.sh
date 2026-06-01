#!/bin/bash

# Check if a commit message was provided
if [ -z "$1" ]; then
    echo "Usage: ./push.sh \"your commit message\""
    exit 1
fi

# Check if the remote is already set, if not, set it up
if ! git remote | grep -q "origin"; then
    echo "Setting up remote origin..."
    git remote add origin https://github.com/KeRon-asm/Protein-Sequence-Analyzer.git
    echo "Remote set successfully"
fi

# Stage all changes
git add .

# Commit with the provided message
git commit -m "$1"

# Push to master
git push origin master

echo "Done! Changes pushed to Protein-Sequence-Analyzer."
