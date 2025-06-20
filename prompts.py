RECOPILAR_INFORMACION ="""
Recibe una consulta en lenguaje natural sobre información genética y extrae los siguientes datos específicos:
- gen1: Nombre del primer gen mencionado
- gen2: Nombre del segundo gen mencionado 
- especie: Especie mencionada (por ejemplo, Homo sapiens)

Ejemplo de entrada:
"Necesito comparar las secuencias de los genes BRCA1 y TP53 en humanos."

Ejemplo de salida extraída:
gen1 = "BRCA1"
gen2 = "TP53"
especie = "Homo sapiens"

PARA TENER EN CUENTA:
 -En caso de que el usuario no suministre la información debes solicitar lo que falta.

Una vez, recopilada la información usa la herramienta descargar_y_fusionar_fastas pasándole como argumentos: gen1, gen2, y especie y
Devolver al usuario la información sobre la secuencia fusionada (longitudes, secuencia total).

"""