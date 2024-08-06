if [ -e env/bin/activate ]; then
	. env/bin/activate
else
	echo "[ERROR]: Cannot find $(env/bin/activate). You should source this script directly from the directory that contains it"
	echo "         Alternatively, you may not have built the venv for this package. Do so with:"
	echo "  $ cd env && python3 -m venv && . bin/activate && pip install -r requirements.txt"
fi
