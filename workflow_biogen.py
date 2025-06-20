from agents import Agent, Runner, FileSearchTool, RawResponsesStreamEvent, Runner, TResponseInputItem, trace, enable_verbose_stdout_logging
from openai.types.responses import ResponseContentPartDoneEvent, ResponseTextDeltaEvent
import asyncio
from fastapi import FastAPI, UploadFile, File
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
from Bio.PDB import PDBList
from pathlib import Path
import os # Asegúrate de que os está importado
from fastapi.middleware.cors import CORSMiddleware
from groq import Groq
from uuid import uuid4



load_dotenv()
openai.api_key = os.getenv("OPENAI_API_KEY")
os.environ["GROQ_API_KEY"] = "GROQ_API_KEY=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"
client = Groq(api_key="GROQ_API_KEY=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")


app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["https://agentomics.vercel.app","http://localhost:3000"], 
    allow_credentials=True,
    allow_methods=["*"],  
    allow_headers=["*"],
)

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
@function_tool
def buscar_pdb_por_nombres_fasta(nombres_fasta: str) -> str:
    """
    Dado un string con uno o varios nombres de archivos de secuencias en formato FASTA (por ejemplo, 'BRCA1.fasta,TP53.fasta'),
    identifica el gen correspondiente a cada archivo y busca en internet los archivos PDB (Protein Data Bank)
    asociados a esos genes.

    Parámetros:
        nombres_fasta (str): Nombres de los archivos FASTA separados por comas. Ejemplo: 'BRCA1.fasta,TP53.fasta'

    Retorna:
        str: Una lista (en string) con los nombres de los archivos PDB a descargar asociados a cada nombre de fasta.
        Si no se encuentra un PDB para algún archivo, lo indica en la respuesta.

    Ejemplo:
        buscar_pdb_por_nombres_fasta('BRCA1.fasta,TP53.fasta')
        'BRCA1: 1JNX, TP53: 2AC0'
        Retorna: "1JNX", "2AC0"
    """
    from openai import OpenAI
    openai_client = OpenAI()        
    response = openai_client.responses.create(
        model="gpt-4.1",
        tools=[{
            "type": "web_search_preview",
            "user_location": {
                "type": "approximate",
                "region": "Latin America",
            }            
                
        }],

        input= nombres_fasta,
    )
    return response.output_text

def descargar_estructuras_pdb(ids, carpeta_destino="pdbs"):
    pdbl = PDBList()
    Path(carpeta_destino).mkdir(exist_ok=True)
    archivos = []
    for pdb_id in ids:
        archivo = pdbl.retrieve_pdb_file(pdb_id, pdir=carpeta_destino, file_format='pdb')
        archivos.append(Path(archivo).resolve())
    return archivos

def generar_link_ngl(pdb_id):
    """Genera un link de visualización en NGL Viewer dado un ID de proteína PDB."""
    pdb_id = pdb_id.lower()
    return f"https://nglviewer.org/ngl/?script=load%20https://files.rcsb.org/download/{pdb_id.upper()}.pdb"

def generar_visualizaciones(ids):
    """Descarga y genera links visuales para cada proteína."""
    descargar_estructuras_pdb(ids)
    links = {pdb_id: generar_link_ngl(pdb_id) for pdb_id in ids}
    return links

@function_tool
def descargar_pdbs(ids: str) -> str:
    """
    Tool que recibe un ID de proteína (o múltiples separados por comas) y retorna
    los enlaces de visualización en NGL Viewer para las estructuras 3D.

    Args:
        ids (str): IDs PDB separados por comas (ej. "1TUP,4HHB")

    Returns:
        str: Enlaces de visualización en formato texto.
    """
    list_ids = [id_.strip().upper() for id_ in ids.split(",")]
    enlaces = generar_visualizaciones(list_ids)

    resultado = "\n".join([f"{id_pdb}: {url}" for id_pdb, url in enlaces.items()])
    return resultado 



generar_pdbs = Agent(
    name="generar_pdbs",
    handoff_description="Agente encargado de la busqueda y descarga de pdbs",
    instructions= """
        Eres un agente especializado en buscar proteinas para esto debes usar la herramienta de: 'buscar_pdb_por_nombre_fasta',
        seguidamente usa la herramienta de 'descargar_pbs' para realizar la descargar y generar el link de visualización de la estrcutura 3D.
        Genera una salida limpia para el usuario con una descripcion y el link de visualizacion.

        ### Ejemplo de Salida Esperada

        ---

        **Estructura de Proteína Encontrada**

        *   **Proteína:** Hemoglobina (Humana)
        *   **Identificador PDB:** `2HHB`
        *   **Descripción:** Deoxyhaemoglobin (haemoglobin in the T state). Esta es la forma de la hemoglobina que no está unida al oxígeno.
        *   **Visualización 3D:** https://molstar.org/viewer/?pdb=2HHB

        ---

        Espero que esta información
    """,
    tools=[buscar_pdb_por_nombres_fasta, descargar_pdbs]
)


@function_tool
def search_in_knowledge(input_question: str) -> str:
    """
    Busca la pregunta ingresada en la base de conocimiento (archivos PDF, etc.).

        Argumentos:
            input_text (str): El texto de entrada con la pregunta para consultar al modelo.

        Retorna:
            str: El texto de salida generado por el modelo.
    """

    from openai import OpenAI
    client = OpenAI()

    response = client.responses.create(
        model="gpt-4o",
        input=input_question,
        tools=[{
            "type": "file_search",
            "vector_store_ids": ["vs_6854c883ffb88191967dfaff2848b3e4"]
        }]
    )

    return response.output_text

# AnswerQuestionsAgent: answers questions about the hotel using available knowledge.
agentic_rag = Agent(
    name="agentic_rag",
    handoff_description='''
Eres un Experto en Genética y Biología Molecular. Tu función es actuar como la principal fuente de conocimiento del sistema, respondiendo a las preguntas del usuario con rigor científico y claridad.

### Tu Protocolo de Respuesta:
1.  **Prioriza la Base de Conocimiento:** Al recibir una pregunta, utiliza la herramienta `search_in_knowledge` como tu primera acción.
 Las respuestas basadas en los documentos proporcionados son la fuente de verdad principal.
2.  **Integra Datos del Contexto:** Si la conversación contiene análisis previos (como estadísticas o resultados de fusión), intégralos en tu respuesta para dar un contexto completo.
3.  **Comunica con Claridad:**
    - Traduce conceptos complejos a un lenguaje accesible, sin perder la precisión científica.
    - Cuando cites información de la base de conocimiento, si es posible, menciónalo sutilmente. 
        Ejemplo: "Según la documentación de referencia, el gen BRCA1 está implicado en..."
    - Si la herramienta `search_in_knowledge` no devuelve información relevante, indícalo y responde 
    basándote en tu conocimiento general como modelo entrenado.
             Ejemplo: "No encontré detalles específicos sobre este punto en mis documentos, pero en el campo de la genética, generalmente se entiende que..."
    ''',
    tools=[search_in_knowledge]
)

#agente enrutador
agente_enrutador_genetico = Agent(
    name="agente_enrutador_genetico",
    handoff_description= "Agente encargado de enrutar solicitudes sobre extracción, análisis, consulta de PDBs y generación de reportes PDF.",
    instructions="""
    Eres el Coordinador del Flujo de Trabajo Bioinformático. Tu rol es recibir la solicitud del usuario, entender su intención y delegar la tarea al agente especialista adecuado. Actúas como el punto de contacto principal.

### Reglas de Enrutamiento:
- **Extracción y Fusión:** Si la solicitud implica obtener datos de genes (ej. "compara BRCA1 y TP53"), delega a `extraer_datos_geneticos`.
- **Análisis Estadístico:** Si se pide analizar una secuencia ya creada (ej. "dame las estadísticas de la fusión ","dame mas informacion de un gen"), delega a `generar_estadisticas`.
- **Estructuras 3D:** Si la consulta es sobre estructuras de proteínas o visualizaciones (ej. "muéstrame la estructura del TP53"), delega a `generar_pdbs`.
- **Preguntas Generales:** Para preguntas conceptuales o de conocimiento (ej. "¿Qué función tiene el gen FOXO3?"), delega a `agentic_rag`.

### Cómo Interactuar:
- **Sé un intermediario, no un simple repetidor:** No devuelvas simplemente la salida del otro agente. Preséntala de forma contextualizada.
    - *Ejemplo:* En lugar de solo mostrar el resultado, di: "El agente de análisis ha completado su tarea. Aquí tienes el resumen:"
- **Gestiona la Incertidumbre:** Si la intención del usuario no es clara, haz una pregunta para clarificar antes de delegar. Ejemplo: "Entiendo que quieres información sobre el gen TP53. ¿Te gustaría obtener su secuencia, un análisis estadístico o su estructura 3D?"
- **Mantén el Enfoque:** Si la solicitud está fuera del ámbito de la bioinformática, informa al usuario de manera cortés. Ejemplo: "Mi especialidad es el análisis genético. No puedo ayudarte con ese tema."
    """,
    handoffs=[extraer_datos_geneticos, generar_estadisticas, agentic_rag, generar_pdbs]
)


#persistencia con redis


class PreguntaRequest(BaseModel):
    pregunta: str

class PreguntaResponse(BaseModel):
    respuesta: str


class TranscripcionResponse(BaseModel):
    text: str
import tempfile

@app.post("/transcribir", response_model=TranscripcionResponse)
async def transcribir(file: UploadFile = File(...)):
    try:
        temp_dir = tempfile.gettempdir()
        audio_id = str(uuid4()) + ".webm"  # Cambia extensión si lo necesitas
        audio_path = os.path.join(temp_dir, audio_id)

        # Guardar archivo temporalmente
        with open(audio_path, "wb") as f:
            f.write(await file.read())
        print("[DEBUG] Audio file saved successfully.")

        # Transcribir con Whisper
        from openai import OpenAI
        client = OpenAI()
        print("[DEBUG] Attempting to transcribe audio with Groq...")
        with open(audio_path, "rb") as audio_file:
            response = client.audio.transcriptions.create(
                model="whisper-large-v3",
                file=audio_file,
                response_format="text",
                language="es",
            )
            transcribed_text = response.strip()

        return TranscripcionResponse(text=transcribed_text)
    except Exception as e:
        print(f"[!] Error:", e)
        return JSONResponse(content={"error": str(e)}, status_code=500)
    

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