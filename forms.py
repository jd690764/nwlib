from django import forms
from numpy import unique
import os
from django.contrib.admin.filters import AllValuesFieldListFilter


class lookupForm(forms.Form):

    #symbol     = forms.TextInput(attrs={'placeholder': 'KRAS', 'label': 'Symbol:'})
    #symbol     = forms.CharField(required=True)
    symbol     = forms.CharField(label    = 'Symbol',
                                 required = True,
                                 widget   = forms.TextInput( attrs={'placeholder': 'KRAS'}))
    org        = forms.ChoiceField(choices  = (),
                                   required = True )
    bait       = forms.MultipleChoiceField( widget  = forms.SelectMultiple( attrs = {'size':'20'}),
                                            choices = ())


    def __init__(self, *args, **kwargs):
        super(lookupForm, self).__init__(*args, **kwargs)

        organisms = (  
            ('9606', 'Human'),
            ('10090', 'Mouse'),
            ('all', 'all'),
        )

        ifilenames = [fn for fn in os.listdir('/mnt/msrepo/ifiles') if fn[-2:] == '.i' and not 'bioplex' in fn and fn[0] != '.' ]
        bs         = sorted(unique( [ fn.split('_')[0].upper() for fn in ifilenames ]))
        bs.insert(0, 'all')
        baits      = list()
        for b in bs:
            baits.append((b, b))
        baits = tuple(baits)

        self.fields['org'].choices = organisms
        self.fields['bait'].choices = baits


class DropdownFilter(AllValuesFieldListFilter):
    template = 'admin/dropdown_filter.html'
