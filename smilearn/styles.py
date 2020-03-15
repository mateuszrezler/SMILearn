from IPython.display import display, HTML

BR = '<br/>'


def arrow_down():
    return f'<br/><big>&#x2B07;<big><br/>'


def b(text):
    return f'<b>{text}</b>'


def code(text):
    return f'<code>{text}</code>'


def fig(counter, caption, br=(1, 2)):
    display(HTML(
        '<br/>'*br[0]+f'<br/><b>Fig. {counter}</b> {caption}'+'<br/>'*br[1]))
    counter += 1
    return counter


def hide_atoms(atoms):
    return {atom: [.75, .75, .75] for atom in atoms}


def molimg(mol, size=(320, 320), **mol_to_image_kwargs):
    display(MolToImage(mol, size, **mol_to_image_kwargs))


def i(text):
    return f'<i>{text}</i>'


def p(*texts):
    display(HTML(''.join(texts)))


def tab(counter, caption, br=(1, 2)):
    display(HTML(
        '<br/>'*br[0]+f'<br/><b>Tab. {counter}</b> {caption}'+'<br/>'*br[1]))
    counter += 1
    return counter

