
PROJNAME = "nademo"

import os, shutil

from seamless.highlevel import Context, Cell, Transformer

ctx = None
webctx = None
save = None

async def load():
    from seamless.metalevel.bind_status_graph import bind_status_graph_async
    import json

    global ctx, webctx, save
    graph = json.load(open("graph/" + PROJNAME + ".seamless"))
    for f in (
        "web/index.html", "web/index.js",
        "web/index-CONFLICT.html", "web/index-CONFLICT.js",
        "web/webform.json", "web/webform-CONFLICT.txt"
    ):
        if os.path.exists(f):
            dest = f + "-BAK"
            if os.path.exists(dest):
                os.remove(dest)
            shutil.move(f, dest)
    ctx = Context()
    ctx.load_vault("vault")
    ctx.set_graph(graph, mounts=True, shares=True)
    await ctx.translation(force=True)

    status_graph = json.load(open("graph/" + PROJNAME + "-webctx.seamless"))

    webctx = await bind_status_graph_async(
        ctx, status_graph,
        mounts=True,
        shares=True
    )
    def save():
        import os, itertools, shutil

        def backup(filename):
            if not os.path.exists(filename):
                return filename
            for n in itertools.count():
                n2 = n if n else ""
                new_filename = "{}.bak{}".format(filename, n2)
                if not os.path.exists(new_filename):
                    break
            shutil.move(filename, new_filename)
            return filename

        ctx.save_graph(backup("graph/" + PROJNAME + ".seamless"))
        webctx.save_graph(backup("graph/" + PROJNAME + "-monitoring.seamless"))
        ctx.save_vault("vault")
        webctx.save_vault("vault")

    print("""Project loaded.

    Main context is "ctx"
    Web/status context is "webctx"

    Open http://localhost:<REST server port> to see the web page
    Open http://localhost:<REST server port>/status/status.html to see the status

    Run save() to save the project
    """)
