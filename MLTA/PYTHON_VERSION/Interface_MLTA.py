import tkinter as tk
from tkinter import *
from tkinter import simpledialog, messagebox, Tk, Frame, Label, LEFT, RIGHT, Button, Entry,font
import customtkinter as cttk
from PIL import ImageTk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import json
import os
import re
import time

class Main_wind_pos:
    ## Abre a tela inicial e estabelece algumas informações a serem usadas no meio do código
    def __init__(self, root_1):
        self.root_1 = root_1
        self.root_1.title("MLTA Interface")
        self.fontePadrao = ("Adobe Myungjo Std M", "11")
        self.fonteBotoes = ("Arial", "10")
        cttk.set_appearance_mode("Dark")
        largura_tela_1 = self.root_1.winfo_screenwidth()
        altura_tela_1 = self.root_1.winfo_screenheight()
        self.largura_tela = largura_tela_1//4
        self.altura_tela = altura_tela_1//3
        self.root_1.geometry(f'1000x500+{self.largura_tela}+{self.altura_tela}')
        
    def interface(self):
        ## Aqui deixa já aberto as seguintes funções
        self.define_titulo()
        self.collect_par()
        self.jason_quantities = []
    def define_titulo(self):
    	## Define o título na tela inicial
        self.primeiroContainer = cttk.CTkFrame(self.root_1, corner_radius=15, fg_color="#2B2B2B")
        self.primeiroContainer.pack(pady=10, padx=25, fill="x")
        texto1= """Enter with the following parameters to run the simulation"""
        self.titulo = cttk.CTkLabel(

        self.primeiroContainer,
        text=texto1,
        font=("Ubuntu", 15, "bold"),
        text_color="white",  # Melhor contraste no fundo escuro
        wraplength=498,  # Mantém a formatação organizada 
        justify="center"
        )
        self.titulo.pack(pady=8, padx=15)
    def collect_par(self):
        self.mainconteiner_root1 = cttk.CTkFrame(root_1)
        self.mainconteiner_root1.pack(pady=10, padx=20)

        
        self.titulo_maincont = cttk.CTkLabel(self.mainconteiner_root1, 
                           text="""Parameters""",
                           justify="center")
        self.titulo_maincont.pack(pady=2, padx=2)
        
        ## Fração volumétrica de partículas
        self.segundoContainer = cttk.CTkFrame(self.mainconteiner_root1)
        self.segundoContainer.pack(pady=10, padx=20,fill="both")
        
        self.volfracLabel = cttk.CTkLabel(self.segundoContainer, 
                              text="Volume fraction of particles (%): ")
        self.volfracLabel.pack(side="left")
        self.volfrac = cttk.CTkEntry(self.segundoContainer, 
                         width=300, 
                         placeholder_text="3.5")
        self.volfrac.pack(side="right")
        
        
        ## Intensidade do campo magnético
        self.terceiroContainer = cttk.CTkFrame(self.mainconteiner_root1)
        self.terceiroContainer.pack(pady=10, padx=20)
        
        self.intfieldLabel = cttk.CTkLabel(self.terceiroContainer, 
                              text="Amplitude of the applied field (A/m) : ")
        self.intfieldLabel.pack(side="left")
        self.intfield = cttk.CTkEntry(self.terceiroContainer, 
                         width=300, 
                         placeholder_text="2.0e+03")
        self.intfield.pack(side="right")

        ## Frequência do campo magnético
        self.quartoContainer = cttk.CTkFrame(self.mainconteiner_root1)
        self.quartoContainer.pack(pady=10, padx=20)
        
        self.magfreqLabel = cttk.CTkLabel(self.quartoContainer, 
                              text="Frequency of the applied field (Hz) : ")
        self.magfreqLabel.pack(side="left")
        self.magfield = cttk.CTkEntry(self.quartoContainer, 
                         width=300, 
                         placeholder_text="2.0e+05")
        self.magfield.pack(side="right")

        ## raio da partícula magnética
        self.quintoContainer = cttk.CTkFrame(self.mainconteiner_root1)
        self.quintoContainer.pack(pady=10, padx=20)
        
        self.radLabel = cttk.CTkLabel(self.quintoContainer, 
                              text="Nominal radius of the particles (m) : ")
        self.radLabel.pack(side="left")
        self.rad = cttk.CTkEntry(self.quintoContainer, 
                         width=300, 
                         placeholder_text="6.0e-09")
        self.rad.pack(side="right")
        
        ## Resultados
        
        self.resultContainer = cttk.CTkFrame(self.root_1)
        self.resultContainer.pack(pady=10, padx=20)
        
        self.temp_label = cttk.CTkLabel(self.resultContainer, text="Temperatura Prevista: ")
        self.temp_label.pack()

        self.time_label = cttk.CTkLabel(self.resultContainer, text="Tempo Previsto: ")
        self.time_label.pack()
        
        #botão de gerar jason
        self.autenticar = cttk.CTkButton(self.mainconteiner_root1, text="Ok",
                           width=150,  # Ajustado para largura em pixels
                           height=40,command=self.gera_jason)
        self.autenticar.pack(side=BOTTOM,pady=5,padx=5,fill="none", expand=False)
        
    def gera_jason(self):
        os.system("rm AIFINAL.py")
        os.system("cp AIFINAL_orig.py AIFINAL.py")
        import json
        inputDict_parameter = {}
        inputDict_parameter["volume"] = self.volfrac.get()
        inputDict_parameter["frequency"] = self.magfield.get()
        inputDict_parameter["magnetic_amplitude"] = self.intfield.get()
        inputDict_parameter["radius"] = self.rad.get()
        
        ## Salva no jason
        json_string = json.dumps(inputDict_parameter, indent=4)
        self.outJson="inputDict_parameters.json"
        self.jason_quantities.append(self.outJson)
        
        with open(self.outJson,"w") as f:
            f.write(json_string)
        f.close()

        ## Monta o botão de confirmação e fecha as telas
        self.confirmationmesh = cttk.CTkToplevel(root_1)
            
        self.confContainermesh = cttk.CTkFrame(self.confirmationmesh)
        self.confContainermesh.pack(pady=10, padx=20)
        
        # Título do container
        self.confirmationmesh.title("Confirmation")
        self.confirmationmeshLabel = cttk.CTkLabel(self.confContainermesh, 
                              text="Got it! Please wait a minute or two for ANN's calculations.")
        self.confirmationmeshLabel.pack(side="left")
        
        ## Botão de confirmação
        self.confirmationmeshbutton = cttk.CTkButton(self.confirmationmesh, text="Ok",
                           width=150,  # Ajustado para largura em pixels
                           height=40,command=self.window_destroymesh)
        self.confirmationmeshbutton.pack(side=BOTTOM,padx=5)
        
    def window_destroymesh(self):
        self.confirmationmesh.destroy()
        self.call_file()
    def call_file(self):
        import os
        json_file_path1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputDict_parameters.json")
        with open(json_file_path1,'r') as f:
            data = json.load(f)
        inputDict = self.generate_dictionary_2(data)
        self.edita_file(inputDict)
        self.run_script()
    def generate_dictionary_2(self,data,dir=None):
        if dir is None:
        	dir = os.path.dirname(os.path.abspath("AIFINAL.py"))
        	dict1 = {
        	os.path.join(dir, "AIFINAL.py"):
        	[
        	{"volume":{"exp":"{phi}","value":data["volume"]}},
        	{"frequency":{"exp":"{freq}","value":data["frequency"]}},
        	{"magnetic_amplitude":{"exp":"{H_0}","value":data["magnetic_amplitude"]}},
        	{"radius":{"exp":"{rad}","value":data["radius"]}}
        	]
        }
        return dict1
        
        
    def edita_file(self,dict1,dir=None):
        for file in dict1.keys():
        	with open(file,"r") as f:
            		entrie=f.read()
        	f.close()
        	dfile = dict1[file]
        	for entrie_line in dfile:
            		for key in entrie_line:
                    		exp = entrie_line[key]["exp"]
                    		value = entrie_line[key]["value"]
                    		entrie = re.sub(rf'{exp}', str(value), entrie)
        	with open(file,"w") as f:
            		f.write(entrie)
        	f.close()
    def run_script(self):
        os.system("python3 AIFINAL.py")
        json_file_path = "resultado.json"
        while not os.path.exists(json_file_path):
                # Aguarda 1 segundo antes de tentar novamente
                time.sleep(1)
    
        # Quando o arquivo JSON for encontrado, leia os dados
        with open(json_file_path, "r") as f:
                data1 = json.load(f)
    
        # Extrai os dados do arquivo JSON
        predicted_temp = data1["predicted_temp"]
        predicted_time = data1["predicted_time"]


        # Atualiza os Labels com os resultados
        self.temp_label.configure(text=f"Predicted temperature: {predicted_temp:.2f} °C")
        self.time_label.configure(text=f"Predicted time: {predicted_time:.2f} s")
    

    
    
root_1 = cttk.CTk()
app = Main_wind_pos(root_1)
zaragui = app # Inicialização

zaragui.interface()

root_1.mainloop()
