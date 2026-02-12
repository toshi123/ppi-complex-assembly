## UniProt alias -> primary mapping

1. Export AC list from sqlite
2. Split into <=100k lines
3. UniProt ID mapping: UniProtKB AC/ID -> UniProtKB, download "TSV (from/to only)"
4. Merge TSVs
5. Convert TSV to alias->primary JSON (local one-off script)
