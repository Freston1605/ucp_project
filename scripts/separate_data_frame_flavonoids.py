import pandas as pd


df = pd.read_csv("data\output\csv\flavonoids_data.csv", delimiter=',')

#Cambiar espacios por "_"
df.columns = df.columns.str.replace(' ', '_')

#Numerar filas que no contienen información
flavonoides_sin_identificador = df[df['PubChem_SID'].isnull()]
cantidad_flavonoides_sin_identificador = flavonoides_sin_identificador.shape[0]
print(f'Número de flavonoides sin identificador en PubChem_SID (Archivo Original): {cantidad_flavonoides_sin_identificador}')
print("Listado de flavonoides sin identificador en PubChem_SID (Archivo Original):")
print(flavonoides_sin_identificador[['ChEBI_Name']])


#Extraer valores SID y CID
df['PubChemSID'] = df['PubChem_SID'].str.extract(r'SID:\s*(\d+)')
df['PubChemCID'] = df['PubChem_SID'].str.extract(r'CID:\s*(\d+)')
print("Columnas creadas: SID y CID")
print(df[['PubChemSID', 'PubChemCID']])
#Eliminar columna antigua
df = df.drop(columns=['PubChem_SID'])

#Corroborar filas sin informacion, tienen que ser las mismas que en un comienzo:
flavonoides_sin_SID = df[df['PubChemSID'].isnull()]
cantidad_flavonoides_sin_SID = flavonoides_sin_SID.shape[0]
print(f'Número de flavonoides sin identificador en Nuevo_PubChem_SID: {cantidad_flavonoides_sin_SID}')
print(flavonoides_sin_SID[['ChEBI_Name']])

flavonoides_sin_CID = df[df['PubChemCID'].isnull()]
cantidad_flavonoides_sin_CID = flavonoides_sin_CID.shape[0]
print(f'Número de flavonoides sin identificador en Nuevo_PubChem_CID: {cantidad_flavonoides_sin_CID}')

print(flavonoides_sin_CID[['ChEBI_Name']])

# Guardar Data Frame
Nuevo_flavonoids_data = "data\output\csv\new_flavonoids_data.csv"
df.to_csv(Nuevo_flavonoids_data, index=False)