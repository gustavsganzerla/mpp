from django.shortcuts import render, redirect
from django.urls import reverse
from django.http import HttpResponse
from .forms import proteinForm, ContactForm

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

from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import json
import requests
import zipfile
from django.conf import settings
import os
import tempfile
import shutil

from django.core.mail import EmailMessage, get_connection
from django.conf import settings

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


def aa_counts(protein):
    allowed_aa = ['A','R','N','D','C','Q','E','G','H','I','K','L','M','F','P','S','T','W','Y','V']
    counts = {el: 0 for el in allowed_aa}

    for el in allowed_aa:
        counts[el] = protein.count(el)

    return counts
    


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
                        

                        aa_percentage = get_amino_acids_percent(protein)
                        aa_percentage_non_zero = {k: round(v, 2) for k, v in aa_percentage.items() if v}
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
                            'aa_counts': aa_counts(protein),
                            'description': description,
                            'sequence': protein,
                            'id':i,
                            'atomic_composition': atomic_composition(protein)
                        })
                        i+=1
                    

                                        
                # Redirect user to another page   
                request.session['output'] = output
                return redirect(reverse('protein_profiler:submitted_prot_char'))
            
            if protein_text:
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

                        aa_percentage = get_amino_acids_percent(protein)
                        aa_percentage_non_zero = {k: round(v, 2) for k, v in aa_percentage.items() if v}
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
                            'aa_counts': aa_counts(protein),
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
                            'Secondary Structure Fraction (helix, turn, sheet)',
                            'Molar Extinction Coefficient (reduced, oxidized)',
                            'Atomic Composition',
                            'Amino Acid Composition'
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
                sequence['charge_at_pH'],
                sequence['secondary_structure_fraction'],
                sequence['molar_extinction_coefficient'],
                sequence['atomic_composition'],
                sequence['get_amino_acids_percent']


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
    output = request.session.get('output', None)
    plt.clf()
    if output:

        if choice == 'length':
            length_data = [item['length'] for item in output]

            plt.hist(length_data)  
            plt.xlabel("Amino acids")
            plt.ylabel("Frequency")

            
            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)


            plot_image = base64.b64encode(buffer.read()).decode()


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

def how_to_use(request):
    return render(request, 'protein_profiler/about.html')

def download(request):
    return render(request, 'protein_profiler/download.html')

def download_zip_file(request):
    zip_file_path = os.path.join(settings.MEDIA_ROOT, 'zip_files', 'mpp_local.zip')

    if not os.path.exists(zip_file_path):
        return HttpResponse('File not found', status=404)
    
    with open(zip_file_path, 'rb') as file:
        response = HttpResponse(file.read(), content_type='application/zip')
        response['Content-Disposition'] = 'attachment; filename="mpp_local.zip"'
        return response
            
def publication(request):
    return render(request, 'protein_profiler/publication.html')

def download_plots(request):
    output = request.session.get('output', None)

    if not output:
        return HttpResponse("No data available to download.", status=404)

    temp_dir = tempfile.mkdtemp()

    #individually plotting each column

    #length
    length_data = [item['length'] for item in output]
    plt.hist(length_data)
    plt.xlabel("Amino acids")
    plt.ylabel("Frequency")
    length_plot_path = os.path.join(temp_dir, 'length_plot.png')
    plt.savefig(length_plot_path)
    plt.close()

    #GRAVY
    gravy = [item['gravy'] for item in output]
    x = [item['id'] for item in output]
    plt.scatter(y = gravy, x = x)
    plt.xlabel("Protein ID")
    plt.ylabel("GRAVY")
    gravy_plot_path = os.path.join(temp_dir, 'gravy_plot.png')
    plt.savefig(gravy_plot_path)
    plt.close()

    #aliphatic index
    aliphatic_index = [item['aliphatic_index'] for item in output]
    x = [item['id'] for item in output]
    plt.scatter(y = aliphatic_index, x = x)
    plt.xlabel("Protein ID")
    plt.ylabel("Aliphatic Index")
    aliphatic_index_plot_path = os.path.join(temp_dir, 'aliphatic_index.png')
    plt.savefig(aliphatic_index_plot_path)
    plt.close()

    #instability index
    instability_index = [item['instability_index'] for item in output]
    x = [item['id'] for item in output]
    plt.scatter(y = instability_index, x = x)
    plt.xlabel("Protein ID")
    plt.ylabel("Instability Index")
    plt.axhline(y = 40, color = 'r', linestyle = '-') 
    instability_index_plot_path = os.path.join(temp_dir, 'instability_index.png')
    plt.savefig(instability_index_plot_path)
    plt.close()

    #stability
    stability = [item['stability'] for item in output]
    stability_counts = Counter(stability)
    labels = stability_counts.keys()
    sizes = stability_counts.values()
    colors = ['green', 'red']
    plt.figure(figsize=(6, 6))
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)
    stability_plot_path = os.path.join(temp_dir, 'stability.png')
    plt.savefig(stability_plot_path)
    plt.close()

    #molecular weight
    molecular_weight = [item['molecular_weight'] for item in output]
    x = [item['id'] for item in output]
    plt.scatter(y = molecular_weight, x = x)
    plt.xlabel("Protein ID")
    plt.ylabel("Molecular Weight")
    molecular_weight_plot_path = os.path.join(temp_dir, 'molecular_weight.png')
    plt.savefig(molecular_weight_plot_path)
    plt.close()

    #aromaticity
    aromaticity = [item['aromaticity'] for item in output]
    x = [item['id'] for item in output]
    plt.scatter(y = aromaticity, x = x)
    plt.xlabel("Protein ID")
    plt.ylabel("Aromaticity")
    aromaticity_plot_path = os.path.join(temp_dir, 'aromaticity.png')
    plt.savefig(aromaticity_plot_path)
    plt.close()

    #isoelectric point
    isoelectric_point = [item['isoelectric_point'] for item in output]
    x = [item['id'] for item in output]
    plt.scatter(y = isoelectric_point, x = x)
    plt.xlabel("Protein ID")
    plt.ylabel("Isoelectric Point")
    isoelectric_point_plot_path = os.path.join(temp_dir, 'isoelectric_point.png')
    plt.savefig(isoelectric_point_plot_path)
    plt.close()

    #charge at ph 7
    charge_at_pH = [item['charge_at_pH'] for item in output]
    x = [item['id'] for item in output]
    plt.scatter(y = charge_at_pH, x = x)
    plt.xlabel("Protein ID")
    plt.ylabel("Charge at pH 7")
    charge_at_pH_plot_path = os.path.join(temp_dir, 'charge_at_pH.png')
    plt.savefig(charge_at_pH_plot_path)
    plt.close()


    # Create a BytesIO buffer to store the zip file contents
    zip_buffer = BytesIO()

    # Create a ZipFile object in memory
    with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
        # Add plot files to the zip file
        zip_file.write(length_plot_path, arcname='length_plot.png')
        zip_file.write(gravy_plot_path, arcname='gravy_plot.png')
        zip_file.write(aliphatic_index_plot_path, arcname='aliphatic_index.png')
        zip_file.write(instability_index_plot_path, arcname='instability_index.png')
        zip_file.write(stability_plot_path, arcname='stability.png')
        zip_file.write(molecular_weight_plot_path, arcname='molecular_weight.png')
        zip_file.write(aromaticity_plot_path, arcname='aromaticity.png')
        zip_file.write(isoelectric_point_plot_path, arcname='isoelectric_point.png')
        zip_file.write(charge_at_pH_plot_path, arcname='charge_at_pH.png')

    # Delete the temporary directory and its contents
    os.remove(length_plot_path)
    os.remove(gravy_plot_path)
    os.remove(aliphatic_index_plot_path)
    os.remove(instability_index_plot_path)
    os.remove(stability_plot_path)
    os.remove(molecular_weight_plot_path)
    os.remove(aromaticity_plot_path)
    os.remove(isoelectric_point_plot_path)
    os.remove(charge_at_pH_plot_path)

    # Seek to the beginning of the BytesIO buffer
    zip_buffer.seek(0)

    # Prepare the HTTP response with the zip file data as an attachment
    response = HttpResponse(zip_buffer.read(), content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename="plots.zip"'
    return response
    
    
def contact(request):
    if request.method=='POST':
        form = ContactForm(request.POST)

        if form.is_valid():
            collected_data = form.cleaned_data
            subject = collected_data.get('subject')
            email = collected_data.get('email')
            message = collected_data.get('message')

            with get_connection(
                host = settings.EMAIL_HOST,
                port = settings.EMAIL_PORT,
                username = settings.EMAIL_HOST_USER,
                password = settings.EMAIL_HOST_PASSWORD,
                use_ssl = settings.EMAIL_USE_SSL
            ) as connection:
                subject = "MPP_"+subject
                email_from = settings.EMAIL_HOST_USER
                recipient_list = ['sganzerlagustavo@gmail.com']
                message = f"{message}\n{email}"

                email = EmailMessage(
                    subject,
                    message,
                    email_from,
                    recipient_list
                )
                email.send()
                return render(request, 'protein_profiler/contact_success.html')
    
    else:
        form = ContactForm()





    return render(request, 'protein_profiler/contact.html', {'form': form})      
