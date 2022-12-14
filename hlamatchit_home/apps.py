from django.apps import AppConfig


class hlamatchitConfig(AppConfig):
    #default_auto_field = 'django.db.models.BigAutoField'
    name = 'hlamatchit_home'

class allele2aa(AppConfig):
    #default_auto_field = 'django.db.models.BigAutoField'
    name = 'allele2aa'

class allele2string(AppConfig):
    #default_auto_field = 'django.db.models.BigAutoField'
    name = 'allele2string'

class ismatch(AppConfig):
    #default_auto_field = 'django.db.models.BigAutoField'
    name = 'ismatch'

class mmcount(AppConfig):
    #default_auto_field = 'django.db.models.BigAutoField'
    name = 'mmcount'


class antigen2aa(AppConfig):
    #default_auto_field = 'django.db.models.BigAutoField'
    name = 'antigen2aa'


class stringcount(AppConfig):
    #default_auto_field = 'django.db.models.BigAutoField'
    name = 'stringcount'