echo "static const char* sources[] ="
echo "{"
for f; do
	echo "\"\\"
	cat $f | while IFS= read -r line; do
		echo "$line\\n\\"
	done
	echo "\","
done
echo "};"
