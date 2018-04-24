## Pangea build

### How to run

* Run all commands from this directory

```
python pangea.prepare.py
```
* Running `pangea.prepare.py` creates the files `prepared/pangea_env.json`, `prepared/pangea_gag.json` and `prepared/pangea_pol.json`.

#### 3. Run build
```
python pangea.process.py --json prepared/pangea_env.json
python pangea.process.py --json prepared/pangea_gag.json
python pangea.process.py --json prepared/pangea_pol.json
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 4. Copy JSONs to auspice
```
cp auspice/pangea_* ../../../auspice/data/
```

#### 5. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
