import pandas as pd


df = pd.read_csv(r"data\output\csv\flavonoids_data.csv", delimiter=',')

#Cambiar espacios por "_"
df.columns = df.columns.str.replace(' ', '_')

#Numerar filas que no contienen información
flavonoides_sin_identificador = df[df['PubChem_SID'].isnull()]
cantidad_flavonoides_sin_identificador = flavonoides_sin_identificador.shape[0]

#Extraer valores SID y CID
df['PubChemSID'] = df['PubChem_SID'].str.extract(r'SID:\s*(\d+)')
df['PubChemCID'] = df['PubChem_SID'].str.extract(r'CID:\s*(\d+)')

#Eliminar columna antigua
df = df.drop(columns=['PubChem_SID'])

#Corroborar filas sin informacion
#Nueva columna "PubChemSID"
flavonoides_sin_SID = df[df['PubChemSID'].isnull()]
cantidad_flavonoides_sin_SID = flavonoides_sin_SID.shape[0]
#Nueva columna "PubChemCID"
flavonoides_sin_CID = df[df['PubChemCID'].isnull()]
cantidad_flavonoides_sin_CID = flavonoides_sin_CID.shape[0]

# Guardar Data Frame
Nuevo_flavonoids_data = r"data\output\csv\new_flavonoids_data.csv"
df.to_csv(Nuevo_flavonoids_data, index=False)