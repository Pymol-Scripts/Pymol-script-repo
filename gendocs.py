import ast
import tempfile
import webbrowser
from collections import defaultdict
from glob import glob
from textwrap import dedent, indent
from jinja2 import Template


FOUND = defaultdict(dict)
NeedRefactor = object()


class Visitor1pass(ast.NodeVisitor):
    def __init__(self, module):
        self.module = module

    def visit_Call(self, call):
        if hasattr(call.func, 'attr') and call.func.attr == "extend":
            if len(call.args) == 2:
                if call.func.attr == "extend" and call.func.value.id == "cmd":
                    funcname = call.args[0].value
                    FOUND[self.module][funcname] = None

class Visitor2Pass(ast.NodeVisitor):
    def __init__(self, module):
        self.module = module

    def visit_FunctionDef(self, functionDef):
        if functionDef.col_offset == 0 and functionDef.name == "__init__":
            FOUND[self.module] = NeedRefactor


class Visitor3pass(ast.NodeVisitor):
    def __init__(self, module):
        self.module = module

    def visit_FunctionDef(self, functiondef):
        if functiondef.name in FOUND[self.module]:
            try:
                if isinstance(functiondef.body[0].value.value, str):
                    docstring = functiondef.body[0].value.value
                    docstring = indent(dedent(docstring), " " * 4)
                    FOUND[self.module][functiondef.name] = docstring
            except Exception as exc:
                pass


for module in [*glob("*.py"), *glob("plugins/*.py")]:
    if module.startswith('_'):
        continue
    with open(module) as module_file:
        src = ''.join(module_file.readlines())
        tree = ast.parse(src, filename=module)

        v = Visitor1pass(module)
        v.visit(tree)

        v = Visitor2Pass(module)
        v.visit(tree)

        if FOUND[module] is not NeedRefactor:
            v = Visitor3pass(module)
            v.visit(tree)
        else:
            del FOUND[module]
            print("%s::NeedRefactor" % module)

tmpl = Template(dedent("""
    {% for module in commands %}
    ## {{ module }}
    {% for command in commands[module] %}
    ### {{ command }}
    ```{{ commands[module][command] }}```
    {% endfor %} 
    {% endfor %}
"""))
docs_path = "command_list.md"
with open(docs_path, "w") as fd:
    fd.write(tmpl.render(commands=FOUND))
