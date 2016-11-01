import os
from fabric.api import env, run, cd

# ===== Usage =====

usage = """

--------
staging       : > fab host_staging deploy:<branch>
prod          : > fab host_prod deploy:<branch>

"""

def help():
    print usage

# ===== hosts ======

def host_staging():
    env.user = 'impd'
    env.hosts = ['cephia.impd.co.za']
    env.code_dir = '/home/cephia'

def host_prod():
    env.user = 'cephia'
    env.hosts = ['cephiadb.incidence-estimation.org']
    env.code_dir = "/home/cephia/shiny-inctools/"

def deploy(branch_name="master"):
    print("   Deploying: ** %s **" % branch_name)
    with cd(env.code_dir):
        run("git reset --hard HEAD")
        run("git fetch origin")
        run("git checkout origin/%s" % branch_name)
        run("git pull origin %s" % branch_name)
        run("sudo service shiny-server restart")

    print("Deployed to: %s" % env.hosts[0])
