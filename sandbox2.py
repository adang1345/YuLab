import io
import json

with io.open("../Disordered Region Data/MobiDB/P35498_disorder.json", "r", encoding='utf8') as f:
    text = f.read()
j_disorder = json.loads(text)
print(j_disorder["consensus"]["long"])