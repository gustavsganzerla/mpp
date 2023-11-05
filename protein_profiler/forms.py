from django import forms

class proteinForm(forms.Form):
    protein_text = forms.CharField(widget=forms.Textarea(attrs={"rows":5, 
                                                                "cols":30,
                                                                "placeholder":">sp|A0A7H0DMZ8|IL18B_MONPV Interleukin-18-binding protein\n MRILFLIAFMYGCVHSYVNAVETKCPNLAIVTSSGEFHCSGCVERMPGFSYMYWLANDMKSDEDTKFIEHLGGGIKEDETVRTTDGGITTLRKVLHVTDTNKFAHYRFTCVLITLDGVSKKNIWLK"}), required=False)
    uploaded_file = forms.FileField(
        label = 'Upload your proteins in .FASTA format', required=False)