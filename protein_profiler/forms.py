from django import forms
from django_recaptcha.fields import ReCaptchaField

class proteinForm(forms.Form):
    protein_text = forms.CharField(widget=forms.Textarea(attrs={"rows":8, 
                                                                "cols":40,
                                                                "placeholder":">sp|A0A7H0DMZ8|IL18B_MONPV Interleukin-18-binding protein\n MRILFLIAFMYGCVHSYVNAVETKCPNLAIVTSSGEFHCSGCVERMPGFSYMYWLANDMKSDEDTKFIEHLGGGIKEDETVRTTDGGITTLRKVLHVTDTNKFAHYRFTCVLITLDGVSKKNIWLK"}), required=False)
    uploaded_file = forms.FileField(
        label = 'Upload your proteins in .FASTA format', required=False)
    


class ContactForm(forms.Form):
    subject = forms.CharField(required=True)
    email = forms.EmailField(required=True)
    message = forms.CharField(widget=forms.Textarea(attrs={"rows":10, 
                                                                "cols":40,
                                                                "placeholder":"Type your message here"}),
                                                                required=True)
    captcha = ReCaptchaField()