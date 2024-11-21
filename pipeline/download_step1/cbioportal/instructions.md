Loop through json files, inside the loop open each with ijson and grab the headers from the first object.

```python
import ijson
import csv
import traceback

error_files = []

# outer loop
for json_file in json_files:
	print(f"Starting file {json_file}")

	# create output file
	tsvfile = open(f"{json_file}_data.tsv", "w")

	# open a file handler 
	file = open(json_file, "r")
	# loop through the json records in that file
	for idx, record in enumerate(ijson.items(file, "items")):
		if idx == 0:
			headers = set(record.keys())
			tsv_writer = csv.DictWriter(tsvfile, delimiter="\t", fieldnames=headers)
			tsv_writer.writeheader()
		try:
			tsv_writer.writerow(record)
		except Exception as e:
			traceback.print_exc()
			print(f"failed to write row {idx} from file {json_file}")
			error_files.append(json_file)
			break

print(f"files with errors: {error_files}")
```
