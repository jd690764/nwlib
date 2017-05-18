from django import forms
from numpy import unique
import os
from django.contrib.admin.filters import AllValuesFieldListFilter
from network.models import Sample

modifs = (('Propionamide', 'Propionamide'),
          ('Oxidation', 'Oxidation'),
          ('Double Oxidation', 'Double Oxidation'),
          ('Deamidated', 'Deamidated'),
          ('Gln->pyro-Glu', 'Gln->pyro-Glu'),
          ('Glu->pyro-Glu', 'Glu->pyro-Glu'),
          ('GlyGly(Ubiquitylation)', 'GlyGly(Ubiquitylation)'),
          ('Carbamoylation', 'Carbamoylation'),
          ('Methyl', 'Methyl'),
          ('Acetyl', 'Acetyl'),
          ('Phosphorylation', 'Phosphorylation'),
          ('Biotinylation', 'Biotinylation'),
          ('APEX2-Biotinylation','APEX2-Biotinylation'),
          ('Unknown', 'Unknown'),)

organisms = (  
    ('9606', 'Human'),
    ('10090', 'Mouse'),
    ('all', 'all'),
)


class lookupForm(forms.Form):

    symbol     = forms.CharField(label    = 'Symbol',
                                 required = True,
                                 widget   = forms.TextInput( attrs={'placeholder': 'KRAS'}))
    org        = forms.ChoiceField(choices  = (),
                                   required = True )
    bait       = forms.MultipleChoiceField( widget  = forms.SelectMultiple( attrs = {'size':'20'}),
                                            choices = ())


    def __init__(self, *args, **kwargs):

        super(lookupForm, self).__init__(*args, **kwargs)

        ifilenames = [fn for fn in os.listdir('/mnt/msrepo/ifiles') if fn[-2:] == '.i' and not 'bioplex' in fn and fn[0] != '.' ]
        bs         = sorted(unique( [ fn.split('_')[0].upper() for fn in ifilenames ]))
        bs.insert(0, 'all')
        baits      = list()
        for b in bs:
            baits.append((b, b))
        baits = tuple(baits)

        self.fields['org'].choices = organisms
        self.fields['bait'].choices = baits

class lookupPtmForm(forms.Form):

    symbol     = forms.CharField(label    = 'Symbol',
                                 required = True,
                                 widget   = forms.TextInput( attrs={'placeholder': 'KRAS'}))
    org        = forms.ChoiceField(choices  = (),
                                   required = True )
    bait       = forms.MultipleChoiceField( widget  = forms.SelectMultiple( attrs = {'size':'20'}),
                                            choices = ())
    modif      = forms.ChoiceField(choices  = (),
                                   required = True )
    expt       = forms.ChoiceField(choices  = (),
                                   required = True )

    def __init__(self, *args, **kwargs):

        super(lookupPtmForm, self).__init__(*args, **kwargs)

        dirs       = os.listdir('/mnt/msrepo/fractionFiles')
        samples    = Sample.objects.filter(ff_folder__in=dirs).values( 'label', 'bait_symbol' )
        bs         = ['all'] + sorted(unique( [ k['bait_symbol'].upper() for k in samples ]))
        baits      = list()
        for b in bs:
            baits.append((b, b))
        baits      = tuple(baits)
        ls         = ['all'] + sorted(unique( [ k['label'].upper() for k in samples ]))
        labels     = list()
        for l in ls:
            labels.append((l, l))
        labels     = tuple(labels)        

        self.fields['org'].choices   = organisms
        self.fields['bait'].choices  = baits
        self.fields['expt'].choices  = labels
        self.fields['modif'].choices = sorted(modifs)
        

class DropdownFilter(AllValuesFieldListFilter):
    template = 'admin/dropdown_filter.html'

