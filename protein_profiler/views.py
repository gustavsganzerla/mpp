from django.shortcuts import render
from django.http import HttpResponse
from .forms import proteinForm

# Create your views here.
def test_view(request):
    return render(request, "protein_profiler/home.html")




def prot_char(request):
    #form = proteinForm()

    if request.method == "POST":
        form = proteinForm(request.POST, request.FILES)
        if form.is_valid():
            collected_data = form.cleaned_data
            protein_text = collected_data.get('protein_text')
            uploaded_file = collected_data.get('uploaded_file')
    else:
        form = proteinForm()

            

    return render(request, 'protein_profiler/prot_char.html', {'form': form})
