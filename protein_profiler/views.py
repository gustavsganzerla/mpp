from django.shortcuts import render, redirect
from django.urls import reverse
from django.http import HttpResponse
from .forms import proteinForm




# Create your views here.
def test_view(request):
    return render(request, "protein_profiler/home.html")




def prot_char(request):

    if request.method == "POST":
        form = proteinForm(request.POST, request.FILES)
        output = []
        if form.is_valid():
            collected_data = form.cleaned_data
            protein_text = collected_data.get('protein_text')
            uploaded_file = collected_data.get('uploaded_file')

            #if uploaded_file:
                #uploaded_file_data = uploaded_file.read().decode('utf-8')
                #request.session['uploaded_file_data'] = uploaded_file_data

                #uploaded_file_io = StringIO(uploaded_file_data)

                #for record in SeqIO.parse(uploaded_file_io, "fasta"):
                    #protein = record.seq
                    #protein = str(protein)

                    #description = record.description

                    #output.append({
                    #    'length':len(protein)
                    #})

                    #request.session['output'] = output
                    #return redirect(reverse('protein_profiler:submitted_prot_char'))
    else:
        form = proteinForm()

            

    return render(request, 'protein_profiler/prot_char.html', {'form': form})


def submitted_prot_char(request):
    output = request.session.get('output', None)

    return render(request, "protein_profiler/submitted_prot_char.html",
                  context={'output':output})