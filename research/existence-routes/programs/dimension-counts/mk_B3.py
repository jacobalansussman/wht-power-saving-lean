import json
sp = json.load(open("specs_B1.json"))[4:] + json.load(open("specs_B2.json"))
for s in sp: s["jac"] = False
json.dump(sp, open("specs_B3.json", "w")); print(len(sp))
