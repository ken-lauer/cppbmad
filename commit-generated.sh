#!/bin/bash
git add $(bash list-generated.sh)
git commit -m "MAINT: regenerate" $(bash list-generated.sh)
