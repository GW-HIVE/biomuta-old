import json

# Assuming the JSON file is named 'sample.json'
with open('/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations/wt_target_2018_pub_wt_target_2018_pub_mutations_wt_target_2018_pub_methylation_all_mutations.json', 'r') as file:
    data = json.load(file)

# Inspect the first object
first_object = data[0]
print("First JSON Object:")
for key, value in first_object.items():
    print(f"Key: {key}, Value: {value}, Type: {type(value)}")

# Inspect the second object
second_object = data[1]
print("\nSecond JSON Object:")
for key, value in second_object.items():
    print(f"Key: {key}, Value: {value}, Type: {type(value)}")