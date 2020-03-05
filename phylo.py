#! /usr/bin/env python
# Theo

# rename internal node with sum leaf name

t = tr('spe_tree', format = 1 )

count=1
for i in t.iter_descendants():
    if not i.is_leaf():
        ll=["".join(["N", str(count)])]
        count += 1
        for j in  i.get_leaves() :
            ll.append(j.name)
        i.name = "_".join(ll)
    else:
        if i.get_distance(t) > 100:
            i.dist = i.dist - (i.get_distance(t) - 100)
        elif i.get_distance(t) < 100:
            i.dist = i.dist + (i.get_distance(t) - 100)

t.name = "World"
t.write(outfile = "spe_rename", format=1, format_root_node=True)
