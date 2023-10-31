from django import forms

class proteinForm(forms.Form):
    protein_text = forms.CharField(widget=forms.Textarea(attrs={"rows":5, "cols":30}), required=False)
    uploaded_file = forms.FileField(
        label = 'Upload your proteins in .FASTA format', required=False)