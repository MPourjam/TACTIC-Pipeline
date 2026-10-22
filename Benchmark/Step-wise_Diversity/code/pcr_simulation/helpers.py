import json


class DictWriter:

    def __init__(self):
        pass

    @staticmethod
    def flatten_dictionary(dic_to_flatten):
        """
        Flattens a dictionary to write in a delimited file.
        """
        flat_dict = {}

        def _flatten(obj, parent_key=''):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    new_key = f"{parent_key}.{k}" if parent_key else k
                    _flatten(v, new_key)
            else:
                flat_dict[parent_key] = obj

        _flatten(dic_to_flatten)
        return flat_dict

    @staticmethod
    def write_delim(dic_to_write, output_file, delimiter) -> str:
        flat_dict = DictWriter.flatten_dictionary(dic_to_write)
        with open(output_file, 'w') as f:
            headers = flat_dict.keys()
            f.write(delimiter.join(headers) + "\n")
            values = [str(flat_dict[h]) for h in headers]
            f.write(delimiter.join(values) + "\n")

        return output_file

    @staticmethod
    def write_tsv(dic_to_write, output_file):
        return DictWriter.write_delim(dic_to_write, output_file, "\t")

    @staticmethod
    def write_csv(dic_to_write, output_file):
        return DictWriter.write_delim(dic_to_write, output_file, ",")
    @staticmethod
    def write_json(dic_to_write, output_file):
        """Writes statistics to a JSON file."""
        with open(output_file, 'w') as f:
            json.dump(dic_to_write, f, indent=4)

        return output_file
