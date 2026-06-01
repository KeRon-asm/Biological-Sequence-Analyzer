#!/bin/bash

# Check if commit msg is provided
if [ -z "$1"]; then
	echo "Usage ./commit_and_push.sh \"my commit message\""
	exit 1 
fi

# Check if git is initialized, if not, setup
if [ ! -d ".git"]; then
	echo "No git repo found, initializing..."

	REPO_NAME=$(basename $PWD)

	git init
	gh repo create "$Biological-Sequence-Analyzer" --public --source=. --remote=origin
	echo "Repo created successfully"
fi

# Stage, commit, push
git add .
git commit -m "$1"
git push origin master

echo "Done! Changes pushed."
