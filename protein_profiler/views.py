from django.shortcuts import render, redirect
from django.urls import reverse
from django.http import HttpResponse
from .forms import proteinForm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from io import StringIO
from io import BytesIO
import re
import csv
import matplotlib
from collections import Counter

matplotlib.use('agg')  # Use the Agg backend

import matplotlib.pyplot as plt
from io import BytesIO
import base64


# Calculations
def gravy(protein):
    gravy_values = {
        'A':  1.8 , 
        'R': -4.5,  
        'N': -3.5,  
        'D': -3.5,  
        'C':  2.5,  
        'Q': -3.5,  
        'E': -3.5,  
        'G': -0.4,  
        'H': -3.2,  
        'I':  4.5,  
        'L':  3.8,  
        'K': -3.9,  
        'M':  1.9,  
        'F':  2.8,  
        'P': -1.6,  
        'S': -0.8,  
        'T': -0.7,  
        'W': -0.9,  
        'Y': -1.3,  
        'V':  4.2  
        }

    sequence_gravy = []
    for aminoacid in protein:
        gravy_value = gravy_values.get(aminoacid)
        sequence_gravy.append(gravy_value)
    return round(sum(sequence_gravy)/len(protein),2)

def aliphatic_index(protein):

    a = 2.9
    b = 3.9

    counts = {
        'A': 0, 
        'V': 0,
        'I': 0, 
        'L': 0, 
    }

    total_amino_acids = len(protein)

    for amino_acid in protein:
        if amino_acid in counts:
            counts[amino_acid] += 1

    mole_percent = {aa: (count / total_amino_acids) * 100 for aa, count in counts.items()}

    aliphatic_index = mole_percent['A'] + a * mole_percent['V'] + b * (mole_percent['I'] + mole_percent['L'])

    return round(aliphatic_index,2)

def instability_index(protein):
    x = ProteinAnalysis(protein)
    return(round(x.instability_index(),2))

def molecular_weight(protein):
    x = ProteinAnalysis(protein)
    return(round(x.molecular_weight(),2))

def aromaticity(protein):
    x = ProteinAnalysis(protein)
    return(round(x.aromaticity(),2))

def isoelectric_point(protein):
    x = ProteinAnalysis(protein)
    return(round(x.isoelectric_point(),2))

def atomic_composition(protein):
    c = 0
    h = 0
    o = 0
    n = 0
    s = 0

    for aa in protein:
        if aa == 'A': #Alanine (Ala)
            c+=3
            h+=7
            o+=2
            n+=1
            s+=0
        if aa == 'R': #Arginine (Arg)
            c+=6
            h+=14
            o+=2
            n+=4
            s+=0
        if aa == 'N':#Asparagine (Asn)
            c+=4
            h+=8
            o+=3
            n+=2
            s+=0
        if aa == 'D':#Aspartic Acid (Asp)
            c+=4
            h+=7
            o+=4
            n+=1
            s+=0
        if aa == 'C':#Cysteine (Cys)
            c+=3
            h+=7
            o+=2
            n+=1
            s+=1
        if aa == 'Q':#Glutamine (Gln)
            c+=5
            h+=10
            o+=3
            n+=2
            s+=0
        if aa == 'E':#Glutamic Acid (Glu)
            c+=5
            h+=9
            o+=4
            n+=1
            s+=0
        if aa == 'G':#Glycine (Gly)
            c+=2
            h+=5
            o+=2
            n+=1
            s+=0
        if aa == 'H':#Histidine (His)
            c+=6
            h+=9
            o+=2
            n+=3
            s+=0
        if aa == 'I':#Isoleucine (Ile)
            c+=6
            h+=13
            o+=2
            n+=1
            s+=0
        if aa == 'L':#Leucine (Leu)
            c+=6
            h+=13
            o+=2
            n+=1
            s+=0
        if aa == 'K':#Lysine (Lys)
            c+=6
            h+=14
            o+=2
            n+=2
            s+=0
        if aa == 'M':#Methionine (Met)
            c+=5
            h+=11
            o+=2
            n+=1
            s+=1
        if aa == 'F':#Phenylalanine (Phe)
            c+=9
            h+=11
            o+=2
            n+=1
            s+=0
        if aa == 'P':#Proline (Pro)
            c+=5
            h+=9
            o+=2
            n+=1
            s+=0
        if aa == 'S':#Serine (Ser
            c+=3
            h+=7
            o+=3
            n+=1
            s+=0
        if aa == 'T':#Threonine (Thr)
            c+=4
            h+=9
            o+=3
            n+=1
            s+=0
        if aa == 'W':#Tryptophan (Trp)
            c+=11
            h+=12
            o+=2
            n+=2
            s+=0
        if aa == 'Y':#Tyrosine (Tyr)
            c+=9
            h+=11
            o+=3
            n+=1
            s+=0
        if aa == 'V':#Valine (Val)
            c+=5
            h+=11
            o+=2
            n+=1
            s+=0


            
    return c,h,o,n,s

def secondary_structure_fraction(protein):
    x = ProteinAnalysis(protein)
    raw_values = x.secondary_structure_fraction()
    rounded_values = tuple(round(value, 2) for value in raw_values)  # Round to 2 decimal places
    return rounded_values

def molar_extinction_coefficient(protein):
    x = ProteinAnalysis(protein)
    return(x.molar_extinction_coefficient())

def charge_at_pH(protein, ph):
    x = ProteinAnalysis(protein)
    return round(x.charge_at_pH(ph),2)

def count_amino_acids(protein):
    x = ProteinAnalysis(protein)
    return(x.count_amino_acids())

def get_amino_acids_percent(protein):
    x = ProteinAnalysis(protein)
    return x.get_amino_acids_percent()




# Create your views here.
def test_view(request):
    return render(request, "protein_profiler/home.html")

def prot_char(request):
    if request.method == 'POST':
        form = proteinForm(request.POST, request.FILES)
        output = []
        stability = ""
        i=1

        if form.is_valid():
            collected_data = form.cleaned_data
            protein_text = collected_data.get('protein_text')
            uploaded_file = collected_data.get('uploaded_file')

            if uploaded_file:
                uploaded_file_data = uploaded_file.read().decode('utf-8')
                
                request.session['uploaded_file_data'] = uploaded_file_data
                
                uploaded_file_io = StringIO(uploaded_file_data)
                
                for record in SeqIO.parse(uploaded_file_io, "fasta"):

                    allowed_amino_acids = set("ARNDCQEGHIJKLMFPSTWYV")

                    if all(aa in allowed_amino_acids for aa in str(record.seq)):

                    
                        protein = record.seq
                        protein = str(protein)
                    
                        
                        description = record.description
                        accession = "NONE"

                        accession = description
                        accession = accession.split(" ")
                        accession = accession[0]
                        
                        if (instability_index(protein)) < 40:
                            stability = "Stable"
                        else:
                            stability = "Unstable"
                        

                        output.append({
                            'accession': accession,
                            'length': len(protein),
                            'gravy': gravy(protein),
                            'aliphatic_index': aliphatic_index(protein),
                            'instability_index': instability_index(protein),
                            'stability':stability,
                            'molecular_weight': molecular_weight(protein),
                            'aromaticity': aromaticity(protein),
                            'isoelectric_point': isoelectric_point(protein),
                            'amino_acid_count': count_amino_acids(protein),
                            'secondary_structure_fraction': secondary_structure_fraction(protein),
                            'molar_extinction_coefficient': molar_extinction_coefficient(protein),
                            'charge_at_pH': charge_at_pH(protein, 7),
                            'get_amino_acids_percent': get_amino_acids_percent(protein),
                            'description': description,
                            'sequence': protein,
                            'id':i,
                            'atomic_composition': atomic_composition(protein)
                        })
                        i+=1
                    

                                        
                # Redirect user to another page   
                request.session['output'] = output
                return redirect(reverse('protein_profiler:submitted_prot_char'))
            
            elif protein_text:
                protein_data = protein_text

                for record in SeqIO.parse(StringIO(protein_data), "fasta"):
                    allowed_amino_acids = set("ARNDCQEGHIJKLMFPSTWYV")

                    if all(aa in allowed_amino_acids for aa in str(record.seq)):

                    
                        protein = record.seq
                        protein = str(protein)
                    
                        description = record.description
                        accession = "NONE"

                        
                            
                        accession = description
                        accession = accession.split(" ")
                        accession = accession[0]


                        if (instability_index(protein)) < 40:
                            stability = "Stable"
                        else:
                            stability = "Unstable"

                        output.append({
                            'accession': str(accession),
                            'length': len(protein),
                            'gravy': gravy(protein),
                            'aliphatic_index': aliphatic_index(protein),
                            'instability_index': instability_index(protein),
                            'stability':stability,
                            'molecular_weight': molecular_weight(protein),
                            'aromaticity': aromaticity(protein),
                            'isoelectric_point': isoelectric_point(protein),
                            'amino_acid_count': count_amino_acids(protein),
                            'secondary_structure_fraction': secondary_structure_fraction(protein),
                            'molar_extinction_coefficient': molar_extinction_coefficient(protein),
                            'charge_at_pH': charge_at_pH(protein, 7),
                            'get_amino_acids_percent': get_amino_acids_percent(protein),
                            'description': description,
                            'sequence': protein,
                            'id':i,
                            'atomic_composition': atomic_composition(protein)
                        })
                        i+=1
                    

                                        
                # Redirect user to another page   
                request.session['output'] = output
                return redirect(reverse('protein_profiler:submitted_prot_char'))

    else:
        form = proteinForm()
    

    return render(request, "protein_profiler/prot_char.html", context={'form': form})


def submitted_prot_char(request):
    output = request.session.get('output', None)

    if output:
        csv_content = []
        csv_content.append(['Accession',
                            'Sequence',
                            'Description',
                            'Length',
                            'GRAVY',
                            'Aliphatic Index',
                            'Instability Index',
                            'Stability',
                            'Molecular Weight',
                            'Aromaticity',
                            'Isoelectric Point',
                            'Charge at pH 7',
                            ])
        
        for sequence in output:
            csv_content.append([
                sequence['accession'],
                sequence['sequence'],
                sequence['description'],
                sequence['length'],
                sequence['gravy'],
                sequence['aliphatic_index'],
                sequence['instability_index'],
                sequence['stability'],
                sequence['molecular_weight'],
                sequence['aromaticity'],
                sequence['isoelectric_point'],
                sequence['charge_at_pH']

            ])
        request.session['csv_content'] = csv_content

        return render(request, "protein_profiler/submitted_prot_char.html", context={'output': output,
                                                                   'csv_content':csv_content})
    else:
        return redirect('protein_profiler:prot_char')
    

def download_csv(request):
    csv_content = request.session.get('csv_content', None)

    if csv_content:
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="output.csv"'

        csv_writer = csv.writer(response)
        for row in csv_content:
            csv_writer.writerow(row)

        return response
    else:
        # Handle the case where there is no CSV content
        return HttpResponse("No CSV content available.")
    
def view_protein(request, id):
    output = request.session.get('output', None)
    print(id)
    protein = []

    if output:
        for sequence in output:
            
            if str(id) == str(sequence['id']):
                
                protein.append({
                    'accession':sequence['accession'],
                    'description':sequence['description'],
                    'length':sequence['length'],
                    'sequence':sequence['sequence'],
                    'gravy':sequence['gravy'],
                    'aliphatic_index':sequence['aliphatic_index'],
                    'instability_index':sequence['instability_index'],
                    'stability':sequence['stability'],
                    'molecular_weight':sequence['molecular_weight'],
                    'aromaticity':sequence['aromaticity'],
                    'isoelectric_point':sequence['isoelectric_point'],
                    'charge_at_pH':sequence['charge_at_pH'],
                    'get_amino_acids_percent': sequence['get_amino_acids_percent'],
                    'atomic_composition':sequence['atomic_composition']

                })

    return render(request, "protein_profiler/view_protein.html", context={'protein':protein})


def plot(request, choice):
    # Retrieve the 'output' data from the session
    output = request.session.get('output', None)
    plt.clf()
    if output:

        if choice == 'length':
            length_data = [item['length'] for item in output]
            # Create the distribution plot
            plt.hist(length_data)  # Adjust the number of bins as needed
            plt.xlabel("Amino acids")
            plt.ylabel("Frequency")

            # Save the plot as an image
            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            # Encode the plot as base64
            plot_image = base64.b64encode(buffer.read()).decode()

            # Pass the plot image to the template
            context = {
                'plot_image': plot_image,
            }
        if choice == 'gravy':
            gravy = [item['gravy'] for item in output]
            x = [item['id'] for item in output]
            plt.scatter(y = gravy, x = x)
            plt.xlabel("Protein ID")
            plt.ylabel("GRAVY")

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }
        if choice == 'aliphatic_index':
            aliphatic_index = [item['aliphatic_index'] for item in output]
            x = [item['id'] for item in output]
            plt.scatter(y = aliphatic_index, x = x)
            plt.xlabel("Protein ID")
            plt.ylabel("Aliphatic Index")

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }
        if choice == 'instability_index':
            instability_index = [item['instability_index'] for item in output]
            x = [item['id'] for item in output]
            plt.scatter(y = instability_index, x = x)
            plt.xlabel("Protein ID")
            plt.ylabel("Instability Index")
            plt.axhline(y = 40, color = 'r', linestyle = '-') 

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }
        if choice == 'stability':
            stability = [item['stability'] for item in output]

            # Count the occurrences of each stability value
            stability_counts = Counter(stability)

            # Prepare the data for the pie chart
            labels = stability_counts.keys()
            sizes = stability_counts.values()
            print(stability_counts)

            # Define colors for the pie chart
            colors = ['green', 'red']

            # Create the pie chart
            plt.figure(figsize=(6, 6))
            plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }
        if choice == 'molecular_weight':
            molecular_weight = [item['molecular_weight'] for item in output]
            x = [item['id'] for item in output]
            plt.scatter(y = molecular_weight, x = x)
            plt.xlabel("Protein ID")
            plt.ylabel("Molecular Weight")

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }
        if choice == 'aromaticity':
            aromaticity = [item['aromaticity'] for item in output]
            x = [item['id'] for item in output]
            plt.scatter(y = aromaticity, x = x)
            plt.xlabel("Protein ID")
            plt.ylabel("Aromaticity")

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }
        if choice == 'isoelectric_point':
            isoelectric_point = [item['isoelectric_point'] for item in output]
            x = [item['id'] for item in output]
            plt.scatter(y = isoelectric_point, x = x)
            plt.xlabel("Protein ID")
            plt.ylabel("Isoelectric Point")

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }
        if choice == 'charge_at_pH':
            charge_at_pH = [item['charge_at_pH'] for item in output]
            x = [item['id'] for item in output]
            plt.scatter(y = charge_at_pH, x = x)
            plt.xlabel("Protein ID")
            plt.ylabel("Charge at pH 7")

            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)

            plot_image = base64.b64encode(buffer.read()).decode()

            context = {
                'plot_image': plot_image,
            }

        return render(request, "protein_profiler/plot.html", context)


def about(request):
    return render(request, 'protein_profiler/about.html')