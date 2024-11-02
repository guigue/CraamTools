import datetime as dt

#####################################################################
#
# "La Fecha de Pascua", Ciencia Hoy, V6, N35.
# Richard L. Branham Jr.
#
# Uso:
# >>> import Pascua
# >>> m,n=pascua.compute(yyyy) # yyyy é o ano a calcular
#
# Original em C.
# Portado a Python por Guigue
#           2021-03-30 (ainda sem tomar a vacina da COVID)
#
#####################################################################

def compute(y=dt.datetime.now().year):
    g = y % 19 + 1				# Aureal Number
    c = y // 100 + 1
    x = 3 * c // 4 - 12
    z = (8 * c + 5) // 25 - 5
    d = 5 * y // 4 - x - 10
    e = (11 * g + 20 + z - x) % 30;

    # Epacta
    if ((e == 25) & (g > 11)) | (e == 24):
         e += 1

    # Meton Cycle
    n = 44 - e
    if (n < 21):
        n += 30
    n = n + 7 - (d + n) % 7
    m='Março'
    if (n>31):
        n-=31
        m = 'Abril'
    print('\n\nA data da Páscoa para {0:4d} é {2:2d} de {1:4s} \n\n'.format(y,m,n))
    return m,n
