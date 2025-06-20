from agents import Agent, Runner, FileSearchTool, RawResponsesStreamEvent, Runner, TResponseInputItem, trace, enable_verbose_stdout_logging
from openai.types.responses import ResponseContentPartDoneEvent, ResponseTextDeltaEvent
import asyncio
from dotenv import load_dotenv
from typing_extensions import TypedDict, Any
from agents import Agent, FunctionTool, RunContextWrapper, function_tool
import asyncio
from prompts import RECOPILAR_INFORMACION
from funciones_workflow import descargar_fasta, fusionar_fasta, limpiar_secuencia, contenido_gc, buscar_patron, secuencia_complementaria, traducir_a_proteina
import csv
from Bio.SeqUtils import gc_fraction
from typing import List, Optional, Dict, Any
from Bio import SeqIO
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel
import openai

import os # Asegúrate de que os está importado
import atexit # Módulo para ejecutar código al salir

load_dotenv()
openai.api_key = os.getenv("OPENAI_API_KEY")


app = FastAPI()
archivos_generados_en_sesion = set()

#tool del agente de extracción de datos y fusión de fastas
# En workflow_biogen.py
os.environ["AGENTS_DISABLE_TRACING"] = "1"
os.environ["AGENTS_TELEMETRY_DISABLED"] = "1"

@function_tool
def descargar_y_fusionar_fastas(gen1: str, gen2: str, especie: str) -> str:
    """
    Descarga, fusiona y retorna resumen. Los archivos generados se limpiarán al final de la sesión.
    """
    try:
        email = "sergio.isidro72@gmail.com"
        especie_clean = especie.replace(' ', '_')

        # Usamos prefijos descriptivos
        fasta1_path = descargar_fasta(gen1, especie, email, prefix=f"{gen1}_{especie_clean}")
        if fasta1_path:
            archivos_generados_en_sesion.add(fasta1_path) # <-- RASTREAR

        fasta2_path = descargar_fasta(gen2, especie, email, prefix=f"{gen2}_{especie_clean}")
        if fasta2_path:
            archivos_generados_en_sesion.add(fasta2_path) # <-- RASTREAR
        
        if not fasta1_path or not fasta2_path:
            return f"No se pudo descargar alguna de las secuencias: {gen1}, {gen2}."

        archivos_individuales = [fasta1_path, fasta2_path]
        
        fusion_path = f"./fastas/{gen1}_{gen2}_{especie_clean}_fusion.fasta"
        sec_fusionada, longitudes = fusionar_fasta(archivos_individuales, output_fasta=fusion_path)
        archivos_generados_en_sesion.add(fusion_path) # <-- RASTREAR

        return (f"Se han descargado y fusionado correctamente las secuencias para los genes {gen1} y {gen2} de {especie}.\n"
                f"- Longitud de {gen1}: {longitudes[0]} bases.\n"
                f"- Longitud de {gen2}: {longitudes[1]} bases.\n"
                f"- Longitud total fusionada: {len(sec_fusionada)} bases.\n\n"
                f"La secuencia fusionada ha sido guardada temporalmente en: {fusion_path}\n"
                f"<!-- Rutas de archivos para análisis posterior: {', '.join(archivos_individuales)} -->")
    except Exception as e:
        return f"[ERROR] {str(e)}"


#agente de extracción de datos y fusión de fastas
extraer_datos_geneticos = Agent(
    name="extraer_datos_geneticos",
    handoff_description="Agente encargado de la extracción de los datos geneticos suministrados en la conversación",
    instructions=RECOPILAR_INFORMACION,
    tools=[descargar_y_fusionar_fastas]
)


#agente de analizar fusiones (AGENTIC CRAG)


#agente de generar estadisticas (adecuar para que las estadisticas tengan en cuenta los fastas, lo mismo las funciones que se usan)
@function_tool
def analizar_fusion_genes(
    fasta_files: List[str],
    export_csv: bool = True
) -> Dict[str, Any]:
    """
    Analiza la fusión de genes a partir de sus archivos FASTA. Guarda el análisis detallado
    codón a codón en un archivo CSV temporal y retorna un resumen global del análisis.
    El archivo CSV generado se eliminará automáticamente al finalizar la sesión.

    Args:
        fasta_files (list[str]): Lista de rutas a los archivos FASTA que se analizarán.
        export_csv (bool, opcional): Si es True, exporta el análisis a un archivo CSV.

    Returns:
        dict: Un diccionario con el RESUMEN GLOBAL, los patrones encontrados y la ruta al
              archivo CSV con el análisis completo.
    """
    # 1. Leer las secuencias y nombres desde los archivos FASTA proporcionados
    secuencias = []
    nombres = []
    for ff in fasta_files:
        try:
            record = SeqIO.read(ff, "fasta")
            secuencias.append(str(record.seq))
            # Heurística para extraer un nombre de gen limpio del ID del FASTA
            nombre_gen = record.id.split('[')[0].split('_')[0]
            nombres.append(nombre_gen)
        except Exception as e:
            return {"error": f"Error leyendo el archivo {ff}: {e}"}

    # 2. Realizar el análisis bioinformático (lógica que ya tenías)
    fusionada = "".join(secuencias)
    fusion_limpia = limpiar_secuencia(fusionada)
    codones = [fusion_limpia[i:i+3] for i in range(0, len(fusion_limpia), 3) if len(fusion_limpia[i:i+3]) == 3]

    origen_codones = []
    for idx, seq in enumerate(secuencias):
        seq_limpia = limpiar_secuencia(seq)
        n_codones = len(seq_limpia) // 3
        origen_codones.extend([nombres[idx]] * n_codones)
    if len(origen_codones) < len(codones):
        origen_codones.extend([nombres[-1]] * (len(codones) - len(origen_codones)))

    analisis_completo = []
    prev_origen = None
    for i, codon in enumerate(codones):
        origen = origen_codones[i]
        cambio = "CAMBIO" if prev_origen and origen != prev_origen else ""
        gc_cod = gc_fraction(codon) * 100
        analisis_completo.append({
            "Posicion": i + 1, "Codon": codon, "Origen": origen,
            "Cambio_origen": cambio, "GC%": float(f"{gc_cod:.1f}")
        })
        prev_origen = origen

    # 3. Preparar el resumen y los patrones (los datos ligeros para el LLM)
    resumen = {
        "longitud_total": len(fusion_limpia),
        "cantidad_codones": len(codones),
        "primer_codon": codones[0] if codones else None,
        "ultimo_codon": codones[-1] if codones else None,
        "gc_total": float(f"{contenido_gc(fusion_limpia):.2f}"),
        "traduccion_proteina": traducir_a_proteina(fusion_limpia),
        "secuencia_complementaria": secuencia_complementaria(fusion_limpia)
    }

    patrones = {
        "inicio_ATG": [p[0] + 1 for p in buscar_patron(fusion_limpia, r"ATG")],
        "parada_TAA_TAG_TGA": [p[0] + 1 for p in buscar_patron(fusion_limpia, r"TAA|TAG|TGA")]
    }

    # 4. Guardar el análisis detallado en un archivo CSV y rastrearlo para su limpieza
    ruta_csv_guardado = "No exportado."
    if export_csv:
        output_folder = "./analyses"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        
        nombre_base_csv = '_'.join(nombres)
        nombre_csv = os.path.join(output_folder, f"{nombre_base_csv}_analisis.csv")

        with open(nombre_csv, "w", newline="", encoding='utf-8') as csvfile:
            fieldnames = ["Posicion", "Codon", "Origen", "Cambio_origen", "GC%"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(analisis_completo)
        
        # ¡Paso clave! Añadir el archivo a la lista de limpieza
        archivos_generados_en_sesion.add(nombre_csv)
        print(f"Análisis detallado guardado temporalmente en: {nombre_csv}")
        ruta_csv_guardado = nombre_csv

    # 5. Retornar solo el resumen ligero al agente
    return {
        "resumen_global": resumen,
        "patrones": patrones,
        "ruta_csv_detalles": ruta_csv_guardado
    }


generar_estadisticas = Agent(
    name="generar_estadisticas",
    handoff_description="Agente encargado de la generación de estadisticas claves para realizar analisis",
    instructions= """Eres un agente experto en bioinformática. Cuando recibas una lista de secuencias génicas fusionadas, tu tarea es analizar la fusión de estos genes a nivel de codón usando la herramienta analizar_fusion_genes.
    
    El contexto de la conversación contendrá las rutas a los archivos FASTA individuales que fueron descargados.
    Busca en el historial una línea comentada que contenga 'Rutas de archivos para análisis posterior' para encontrar la lista de archivos.
    Extrae esa lista de rutas y pásala como el argumento `fasta_files` a la herramienta `analizar_fusion_genes`.

    No inventes secuencias. Utiliza únicamente los archivos proporcionados. Responde estructurando el resultado según lo que retorne la tool.
    
    Utiliza la tool analizar_fusion_genes para:
    - Determinar el origen de cada codón.
    - Detectar cambios de fragmento (transición entre genes).
    - Calcular el porcentaje GC por codón y global.
    - Traducir la secuencia fusionada a proteína hasta el primer codón de parada.
    - Identificar posiciones de codón de inicio (ATG) y de parada (TAA, TAG, TGA).
    - Generar un resumen global y, si se solicita, exportar los datos a CSV.

    Responde estructurando el resultado según lo que retorne la tool, sin inventar ni omitir información.
    """,
    tools=[analizar_fusion_genes]
)

#agente de consultar pdbs y retornar pdbs (busqueda en la web)

#agente de analizar información de los otros agentes- generar pdf's (AGENTIC RAG)

#agente enrutador
agente_enrutador_genetico = Agent(
    name="agente_enrutador_genetico",
    handoff_description= "Agente encargado de enrutar solicitudes sobre extracción, análisis, consulta de PDBs y generación de reportes PDF.",
    instructions="""
    Eres un agente enrutador bioinformático.  
    - Si el usuario solicita descargar, extraer o fusionar genes, usa 'extraer_datos_geneticos'.
    - Si solicita análisis estadístico de una secuencia fusionada, usa 'generar_estadisticas'.
    - Si solicita estructuras de proteínas, usa 'consultar_pdbs'.
    - Si pide un reporte en PDF con los resultados, usa 'generar_pdf'.

    Debes devolver exactamente lo que retorne el agente especializado seleccionado, sin agregar ni omitir información.
    En caso de que la información no sea clara o no se entienda que se debe realizar debes preguntar al usuario.
    En caso de que la información no tenga relación con bioinformática y el alcance de los agentes debes expresar que no es un tema para tratar.
    """,
    handoffs=[extraer_datos_geneticos, generar_estadisticas]
)


#endpoints y desplegar en render o railway

import inspect

#persistencia con redis


class PreguntaRequest(BaseModel):
    pregunta: str

class PreguntaResponse(BaseModel):
    respuesta: str

@app.post("/preguntar", response_model=PreguntaResponse)
async def preguntar(request: PreguntaRequest):
    try:
        input_items = [{"role": "user", "content": request.pregunta}]
        result = await Runner.run(agente_enrutador_genetico, input_items)
        return PreguntaResponse(respuesta=result.final_output)
    except Exception as e:
        print(f"[!] Error:", e)
        return JSONResponse(content={"error": str(e)}, status_code=500)
    

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("workflow_biogen:app", host="0.0.0.0", port=8000, reload=True)