import json
import genson
import json_schema_descriptions as jsd

current_schema = None
output_file = "output.json"
new_schema = "new_odmts_schema.json"

with open(output_file) as file:
    output = json.load(file)

builder = genson.SchemaBuilder()

# Start with current schema
if current_schema != None:
    with open(current_schema) as file:
        builder.add_schema(json.load(file))

# Include patterns and null types that may not be detected directly
builder.add_schema({'properties': {
    'parameters': {},
    'design': {'properties':
                   {'hubs': {'patternProperties': {r'^.*$': {}}},
                    'legs': {'patternProperties': {r'^.*$': {}}}}
               },
    'routes': {'patternProperties': {r'^.*$': {}}},
    'trip_splitting': {'anyOf': [
        {'type': 'null'},
        {'type': 'object',
         'patternProperties': {r'^.*$':
             {'patternProperties': {r'^.*$': {
                 'anyOf': [{'type': 'null'},
                           {'type': 'object'}]
             }}}}}
    ]}
}})

# Add properties to schema that are in the output file
builder.add_object(output)

# Add descriptions
schema = builder.to_schema()
jsd.add_descriptions(schema)

# Rebuild the schema and write to file
builder = genson.SchemaBuilder();
builder.add_schema(schema);
with open(new_schema, "w") as file:
    print(builder.to_json(indent=4), file=file)

# Validate output file with the schema
# validator = jsonschema.Draft7Validator(builder.to_schema())
# validator.is_valid(output)