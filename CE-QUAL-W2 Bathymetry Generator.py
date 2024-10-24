'''This program is used to generate bathymetry for CE-QUAL-W2 model using topo contours and centerline of the river. Users need to provide topo contours, centerline of the river, and W2 segments.
Topo lines in contour shapefile should not have geometry error such as dangling lines, self-intersecting lines, etc.
The W2 segments should not be interescted by two centerlines. 
Contour shapefile has Elevation attribute which is used to generate bathymetry.
W2 segments and centerline  shapefile must have 'Segment' column.
Use W2 activate segment numbers for centerline and W2 segments shapefile, e.g., the first segment should be 2, the second segment should be 3, etc.
Do not include inacitve cell numbers the program will add inactive cells to the output.
For inquiries or to report bugs, please contact afshin.shabani64@gamil.com
'''

import geopandas as gpd
import fiona
from shapely.geometry import shape, mapping, MultiLineString, Polygon, LineString
from shapely.ops import linemerge, polygonize
from collections import defaultdict
import os 
import math
import pandas as pd 
import numpy as np 
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt

################ input locations
import tkinter as tk
from tkinter import filedialog, messagebox
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
def convert_multilinestring_to_polygon(input_shapefile, output_shapefile):
    # Read the MultiLineString shapefile
    with fiona.open(input_shapefile, 'r') as input:
        # Create a new schema for the output Polygon shapefile
        schema = {
            'geometry': 'Polygon',
            'properties': {**input.schema['properties']},
        }
        
        # Collect polygons for each feature
        polygons = []
        for feature in input:
            geom = shape(feature['geometry'])
            if isinstance(geom, MultiLineString):
                merged = linemerge(geom)
                poly = polygonize(merged)
                for p in poly:
                    polygons.append({
                        'geometry': mapping(p),
                        'properties': feature['properties']
                    })
        
        # Write the Polygon shapefile
        with fiona.open(output_shapefile, 'w', driver='ESRI Shapefile', crs=input.crs, schema=schema) as output:
            for polygon in polygons:
                output.write(polygon)



def browse_file(entry):
    filename = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, filename)

def browse_directory(entry):
    directory = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, directory)


def run_program():
    try:
        # Example of updating the progress bar
        # update_progress_bar(progress, 0)
        update_progress_bar(progress, progress_label, 0)
        Bathymetery_contours = Bathymetery_entry.get()
        W2_segments = gpd.read_file(W2_segments_entry.get())
        Center_line = center_line_entry.get()
        zone = zone_var.get().split('-')[1]
        contour_intervals = int(contour_intervals_entry.get())
        scale = float(scale_entry.get())
        friction=float(friction_entry.get())
        method=str(volume_method_var.get())
        interval = round(contour_intervals * scale, 1)
        name = name_entry.get()

        output = output_entry.get()
        if not os.path.exists(output):
            os.makedirs(output)

        contours = os.path.join(output, 'contours')
        if not os.path.exists(contours):
            os.makedirs(contours)

        contour_polygon = os.path.join(output, 'contours_polygon')
        if not os.path.exists(contour_polygon):
            os.makedirs(contour_polygon)

        centerlines = os.path.join(output, 'center_lines')
        if not os.path.exists(centerlines):
            os.makedirs(centerlines)

        segment_contour = os.path.join(output, 'Segment_contour')
        if not os.path.exists(segment_contour):
            os.makedirs(segment_contour)

        W2_input = os.path.join(output, 'W2_input')
        if not os.path.exists(W2_input):
            os.makedirs(W2_input)

        # update_progress_bar(progress, 20)
        update_progress_bar(progress, progress_label, 20)
        #  step 1 
        bathymetry=gpd.read_file(Bathymetery_contours)
        bathymetry= bathymetry.to_crs(zone)  # Convert to UTM projection (Zone 12N)
        bathymetry.columns=['Elevation','geometry']

        # Dissolve contours in bathymetry using Z_MAX column value
        dissolved_bathymetry = bathymetry.dissolve(by='Elevation')
        dissolved_bathymetry.reset_index(inplace=True)

        #### Create contours from dissolved bathymetry##############
        for Elevation in dissolved_bathymetry['Elevation'].unique():
            sub = dissolved_bathymetry.loc[dissolved_bathymetry['Elevation'] == Elevation]
            sub.to_file(os.path.join(contours,'%s.shp' % str(Elevation)))

    ################## Convert MultiLineString to Polygon################
        # Step 2
        # update_progress_bar(progress, 40)
        for file in os.listdir(contours):
            if file.endswith('.shp'):
                print(file)
                poly_contour=convert_multilinestring_to_polygon(os.path.join(contours,file),os.path.join(contour_polygon,file))

    ######## Read Lake Center Line ########################
        centerline = gpd.read_file(Center_line)
        centerline = centerline.to_crs(zone)  # Convert to UTM projection (Zone 12N)
        W2_segments = W2_segments.to_crs(zone)
        W2_segments = W2_segments.dissolve(by='Segment')
        W2_segments.reset_index(inplace=True)
    # centerline['length'] = (centerline.length)/1000 # Convert length to km


    ############ Split centerline into smaller polylines using distance using W2 segments #################
        split_centerline=gpd.overlay(centerline,W2_segments,how='intersection')
        split_centerline.to_file(os.path.join(centerlines,'center_lines_split.shp'))
        lines_gdf=split_centerline.copy()
        lines_gdf['length'] = lines_gdf.length
    
    #### Radian for each segment ############################
    # lines_gdf['Segment']=[x+2 for x in range(len(lines_gdf))]
        length=lines_gdf.length
        crid=lines_gdf.Segment
        radian=[] ### cross section radian 
        for Id in lines_gdf.Segment.unique():
            sub=lines_gdf.loc[lines_gdf['Segment']==Id]
            bounds=sub.geometry.bounds
            rad=math.atan2((bounds['maxy']-bounds['miny']),bounds['maxx']-bounds['minx'])
            radian.append(rad)

    ########### Calculate Width ############################
        for fl in os.listdir(contour_polygon):
            if fl.endswith('.shp'):
                clip=gpd.clip(W2_segments,gpd.read_file(os.path.join(contour_polygon,fl)))
                if len(clip)>0:
                    clip_reproject = clip.to_crs(zone)  # Convert to UTM projection (Zone 12N)
                    clip_reproject.to_file(os.path.join(segment_contour,fl)) 
        # Step 3
        # update_progress_bar(progress, 60)
        update_progress_bar(progress, progress_label, 60)
        lst_contour=[]
        for fl in os.listdir(segment_contour):
            if fl.endswith('.shp'):
                print(fl)
                contour_width=gpd.read_file(os.path.join(segment_contour,fl))
                #crosw['Length']=crosw.length
                contour_width['Area']=contour_width.area
                contour_width['Elevation']=[float(fl.replace('.shp','')) for x in range (len(contour_width))]
                contour_width=contour_width[['Segment','Area','Elevation']]
                lst_contour.append(contour_width)
        dfcross=pd.concat(lst_contour)
        dfcross.sort_values(['Segment','Elevation'],inplace=True)
        group=dfcross.groupby(['Segment','Elevation']).sum()
        group.reset_index(inplace=True)

        lstwidth=[]
        elvmin=[]
        elvmax=[]
        for IND,cros in enumerate(group.Segment.unique()):
            dfcr=group.loc[group['Segment']==cros]
            dfcr.reset_index(inplace=True)
            dfcr.drop('index',axis=1,inplace=True)
        #### Make sure value of area increases by width 
            for i in range(1, len(dfcr)):
                if dfcr.loc[i, 'Area'] < dfcr.loc[i - 1, 'Area']:
                    dfcr.loc[i, 'Area'] = dfcr.loc[i - 1, 'Area']
            rowwidth=[] ## store width for each elevation 
            if len(dfcr)==1: ## if there is only one elevation in segment 
                rvolume=interval*dfcr['Area'].values[0]
                dfcr['width']=rvolume/(interval*lines_gdf.loc[lines_gdf['Segment']==cros].length.values[0])
                dfcr['Elevation']=dfcr['Elevation']#+contour_intervals
                lstwidth.append(dfcr)
                elvmin.append(min(dfcr['Elevation'])) # append min elevation for each segment
            else:
                for index, rows in enumerate(dfcr.iterrows()):
                    if index==0:
                        rvolume=interval*dfcr['Area'].values[0]
                        rwidth=rvolume/(interval*lines_gdf.loc[lines_gdf['Segment']==cros].length.values[0])
                    else:
                        # if index < len(dfcr) - 1:
                        if method=='Trapezoidal':
                            rvolume=interval*0.5*(dfcr.iloc[index]['Area']+dfcr.iloc[index-1]['Area']) # trapsoidal volume
                            rwidth=rvolume/(interval*lines_gdf.loc[lines_gdf['Segment']==cros].length.values[0])
                        elif method=='Prismoidal':
                            rvolume=(interval/6)*(dfcr.iloc[index]['Area']+4*np.sqrt(dfcr.iloc[index]['Area']*dfcr.iloc[index-1]['Area'])+dfcr.iloc[index-1]['Area']) # Prismidal volume
                            rwidth=rvolume/(interval*lines_gdf.loc[lines_gdf['Segment']==cros].length.values[0])
                    if rwidth<5 and rwidth!=0: ### This is a thershold for width because W2 become unstable using narrow width 
                        rwidth=10
                        
                    rowwidth.append(rwidth)  
                dfcr['width']=rowwidth
                lstwidth.append(dfcr)
                elvmin.append(min(dfcr['Elevation'])) # append min elevation for each segment 
            if IND==0:
                elvmax.append(max(dfcr['Elevation'])) # append max elevation for most upstream segment 
        dfwidth=pd.concat(lstwidth)
        dfwidth.drop('Area',axis=1,inplace=True)
        
        ## Generate layers 
        bth=pd.DataFrame({'Elevation':[x for x in range(int(min(elvmin)-contour_intervals*5),(int(max(elvmax))+contour_intervals),contour_intervals)]}) ## add inactive layer to the bottom, you can change value of contour_intervals*5
        #bth.to_csv(os.path.join(W2_input,'Elevation_inputs.csv'))
        bth.sort_values('Elevation',ascending=False,inplace=True)
        bth.reset_index(inplace=True)
        bth.drop('index',inplace=True,axis=1)
            # Step 4
        # update_progress_bar(progress, 80)
        update_progress_bar(progress, progress_label, 80)
    ###### organize vertical layers for each segments to combined all of them in 2D array 
        lstbth=[]
        size=[]
        for ID in dfwidth.Segment.unique():
            sdfw=dfwidth.loc[dfwidth['Segment']==ID]
            sdfw=sdfw.loc[sdfw['Elevation']<=max(elvmax)] # This will exclude elevation above the top layer or upstream layer
            delv=max(sdfw['Elevation']) # get the maximum available elevation 
            msdf=pd.merge(bth,sdfw,on='Elevation',how='outer')
            msdf.loc[msdf['Elevation']>delv,'width']=max(sdfw['width'])
            msdf.drop(['Segment','Elevation'],axis=1,inplace=True)
            msdf.interpolate('linear',limit=10,inplace=True,limit_direction='backward')
            msdf.fillna(0,inplace=True)
            # msdf.reset_index(inplace=True)
            lstbth.append(msdf)
            size.append(len(msdf))
            # msdf.to_csv(r'C:\Work\Snake-hell-canyon\Temp\wb.csv',sep='\t')

        dfbth=pd.concat(lstbth,axis=1)

        ### add first and last layer 
        NCL=len(dfbth) #
        flayer=pd.DataFrame([0 for x in range(NCL)])
        nlayer=pd.DataFrame([0 for x in range(NCL)])

        boundary_segments=[x+1 for x in range (0,max(W2_segments.Segment)) if x+1 not in W2_segments.Segment.unique()]
        for bnd in boundary_segments:
            dfbth.insert(bnd-1,str(bnd), flayer)

        # dfbth=pd.concat([flayer,dfbth],axis=1)
        dfbth=pd.concat([dfbth,nlayer],axis=1)

        ## add vertical resolution 
        vr=pd.DataFrame([interval for x in range(NCL)])
        dfbth=pd.concat([vr,dfbth],axis=1)


        ## add top layer 
        MCL=len(dfbth.columns)
        top=pd.DataFrame([float(0) for x in range(MCL)]).T.values
        top[0][0] =float(interval)
        ndf=dfbth.values
        ndf=np.concatenate((top,ndf),axis=0)


        Layer=pd.DataFrame(['LAYERH']+[x for x in range(MCL-1)]).T.values
        ndf=np.concatenate((Layer,ndf),axis=0)

        FRIC=pd.DataFrame(['FRIC']+[friction for x in range(MCL-1)]).T.values
        ndf=np.concatenate((FRIC,ndf),axis=0)

        for bnd in boundary_segments:
            radian.insert(bnd-1,0)
        radian.insert(len(radian)+1,0)
        PHIO=pd.DataFrame(['PHIO']+radian).T.values
        ndf=np.concatenate((PHIO,ndf),axis=0)


        ELWS=pd.DataFrame(['ELWS']+[np.nan for x in range(MCL-1)]).T.values
        ndf=np.concatenate((ELWS,ndf),axis=0)

        Length=lines_gdf.length.tolist()
        for bnd in boundary_segments:
            Length.insert(bnd-1,0)
        Length.insert(len(Length)+1,0)

        DLX=pd.DataFrame(['DLX']+Length).T.values #+[0]
        ndf=np.concatenate((DLX,ndf),axis=0)

        SEG=pd.DataFrame(['SEG']+[x+1 for x in range(MCL-1)]).T.values
        ndf=np.concatenate((SEG,ndf),axis=0)

        Head=pd.DataFrame(['$%s'%name]+[np.nan for x in range(MCL-1)]).T.values
        ndf=np.concatenate((Head,ndf),axis=0)

        dndf=pd.DataFrame(ndf)

        elevation=[*[np.nan for x in range(8)],*(np.sort(np.round([x*scale for x in bth['Elevation']],2)))[::-1]]
        dndf.insert(len(dndf.columns), '', elevation)
        dndf.to_csv(os.path.join(W2_input,'%s.csv'%name),sep=',',header=False,index=False)
        ########## QC 
        QC_bath=dndf.iloc[8:,2:-1]
        QC_bath['Elevation']=elevation[8:]
        QC_bath.index=[x for x in range (len(QC_bath))]
        QC_DLX=dndf.iloc[2,:-1]
        QC_DLX=QC_DLX.T

        def separate_branches(segment_list):
            branches = []
            current_branch = []

            for i in range(len(segment_list)):
                if i == 0:
                    current_branch.append(segment_list[i])
                else:
                    if segment_list[i] - segment_list[i - 1] > 1:
                        branches.append(current_branch)
                        current_branch = []
                    current_branch.append(segment_list[i])
            
            if current_branch:
                branches.append(current_branch)
            
            return branches

        segment_list = W2_segments['Segment']
        branches = separate_branches(segment_list)
        flattened_list = [item for sublist in branches for item in sublist]
        depth=contour_intervals*scale
        ###### append branch number to each segment
        branch_number = []
        for i in range(len(branches)):
            sub_branch=branches[i]
            sub_branch=pd.DataFrame({'Segment':sub_branch})
            sub_branch['Branch']=[i+1 for x in range(len(sub_branch))]
            branch_number.append(sub_branch)
        branch_number=pd.concat(branch_number)

        ###### add elevation to the branches dataframe 
        branch_list=[]
        for i in flattened_list:
            sub_QC_bath=QC_bath[i]
            branch_list.append(sub_QC_bath)

        df_branch=pd.concat(branch_list,axis=1)
        df_branch['Elevation']=QC_bath['Elevation']

        ######### calculate area and volume for each branch
        Segment=[]
        Elevation=[]
        Area=[]
        Volume=[]
        for columns in df_branch.columns[:-1]:
            print(columns)
            sub_df_branch=df_branch[[columns,'Elevation']]
            sub_df_branch=sub_df_branch.loc[sub_df_branch[columns]>0]
            if len (sub_df_branch)==1:
                volume=sub_df_branch[columns]*depth*QC_DLX[columns]
                area=volume/depth
                Segment.append(columns)
                Elevation.append(sub_df_branch['Elevation'].values[0])
                Area.append(area.values[0])
                Volume.append(volume.values[0])
            else:
                for rows in range(len(sub_df_branch)-1):
                    volume=sub_df_branch[columns].iloc[rows]*depth*QC_DLX[columns]
                    area=volume/depth
                    Segment.append(columns)
                    Elevation.append(sub_df_branch['Elevation'].iloc[rows])
                    Area.append(area)
                    Volume.append(volume)

        QC_Area_Volume=pd.DataFrame({'Segment':Segment,'Elevation':Elevation,'Area':Area,'Volume':Volume})
        AreaT=QC_Area_Volume.groupby('Elevation')['Area'].sum()   
        VolumeT=QC_Area_Volume.groupby('Elevation')['Volume'].sum()
        cumsum_volume=VolumeT.cumsum()

        merge=pd.merge(AreaT,cumsum_volume,on='Elevation')
        merge.columns=['Area(m^2)','Volume(m^3)']
        if scale==0.3048:
            merge['Elevation(ft)']=merge.index*3.2808399
            merge['Area (acres)']=merge['Area(m^2)']*0.00024710
            merge['Volume (acre-ft)']=merge['Volume(m^3)']*0.000810714
        merge.to_excel(os.path.join(W2_input,'Area_Volume.xlsx'))
        update_progress_bar(progress, progress_label,100)
        messagebox.showinfo("Success", "Program completed successfully!")
        def show_plot():
            plt.figure(figsize=(10, 14))
            for i in (branch_number['Branch'].unique()):
                plt.subplot(len(branch_number['Branch'].unique()), 1, i)
                columns=branch_number.loc[branch_number['Branch']==i]['Segment']
                sub_branch=QC_bath[columns]
                distance=QC_DLX[columns]
                cumsum_distance=distance.cumsum()
                cumsum_distance=[round(x/1000,2) for x in cumsum_distance]
                sub_QC_bath=QC_bath[columns]
                bathymetry=sub_QC_bath.iloc[:,:].astype(float).to_numpy()
                y=QC_bath['Elevation'].to_numpy()
                x,y=np.meshgrid(cumsum_distance,y)
                plt.pcolormesh(x, y, bathymetry, shading='auto', cmap='terrain',alpha=1.0)#alpha=1,edgecolor='white',,edgecolor='grey' gist_heat
                plt.title(f'Branch {i}')
                plt.colorbar(label='Width (m)')
                # plt.xlabel('Distance (m)')
                # plt.ylabel('Elevation (m)')
                plt.subplots_adjust(hspace=0.3) 
                plt.ylabel('Elevation (m)')
                plt.savefig(os.path.join(W2_input,'%s.jpg'%name),dpi=300, bbox_inches='tight')
            plt.xlabel('Distance (km)')
            root = tk.Tk()
            root.title("Plot")
            # custom_icon = tk.PhotoImage(file=os.path.join(os.getcwd(),'Tetra-Tech-Logo.png'))  # Update the file path to the correct location of the custom icon
            # root.iconphoto(False, custom_icon)
            canvas = FigureCanvasTkAgg(plt.gcf(), master=root)
            canvas.draw()
            canvas.get_tk_widget().pack()
            root.mainloop()

        show_plot()
        plt.close()
         # Step 5
        # update_progress_bar(progress, 100)
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")
        update_progress_bar(progress, 0)

def browse_file(entry):
    file_path = filedialog.askopenfilename(filetypes=[("Shapefiles", "*.shp")])
    if file_path:
        entry.delete(0, tk.END)
        entry.insert(0, file_path)

root = tk.Tk()
def update_progress_bar(progress, progress_label, value):
    progress['value'] = value
    progress_label.config(text=f"{value}%")
    root.update_idletasks()

root.title("CE-QUAL-W2 Bathymetry Generator")

# Set the custom icon
custom_icon = tk.PhotoImage(file='Lake.png')  # Change this to the path of your custom icon
root.iconphoto(False, custom_icon)

# Create a frame to hold all the widgets
frame = tk.Frame(root, padx=10, pady=10)
frame.pack(fill=tk.BOTH, expand=True)

# Create a progress bar widget
progress = ttk.Progressbar(root, orient="horizontal", length=300, mode="determinate")
progress.pack(fill=tk.X, padx=(10, 10), pady=(10, 10))

# Create a label to display the progress percentage
progress_label = tk.Label(root, text="0%")
progress_label.pack()

# Configure grid layout to resize with window
frame.columnconfigure(1, weight=1)

tk.Label(frame, text="Topo Contours (.shp):").grid(row=0, column=0, sticky=tk.W, padx=(0, 5), pady=(5, 5))
Bathymetery_entry = tk.Entry(frame)
Bathymetery_entry.grid(row=0, column=1, sticky=tk.EW, pady=(5, 5))
tk.Button(frame, text="Browse", command=lambda: browse_file(Bathymetery_entry)).grid(row=0, column=2, sticky=tk.EW, padx=(5, 0), pady=(5, 5))

tk.Label(frame, text="W2 Segments (.shp):").grid(row=1, column=0, sticky=tk.W, pady=(5, 5))
W2_segments_entry = tk.Entry(frame)
W2_segments_entry.grid(row=1, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))
tk.Button(frame, text="Browse", command=lambda: browse_file(W2_segments_entry)).grid(row=1, column=2, sticky=tk.EW, pady=(5, 5))

tk.Label(frame, text="River/Branch Center Lines (.shp):").grid(row=2, column=0, sticky=tk.W, pady=(5, 5))
center_line_entry = tk.Entry(frame)
center_line_entry.grid(row=2, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))
tk.Button(frame, text="Browse", command=lambda: browse_file(center_line_entry)).grid(row=2, column=2, sticky=tk.EW, pady=(5, 5))

utm_zones = [f"Zone {i}-EPSG:326{str(i).zfill(2)}" for i in range(1, 61)]
zone_var = tk.StringVar(root)
zone_var.set(utm_zones[11])  # set default value to Zone 12 (EPSG:32612)

tk.Label(frame, text="UTM Zone:").grid(row=3, column=0, sticky=tk.W, pady=(5, 5))
zone_dropdown = tk.OptionMenu(frame, zone_var, *utm_zones)
zone_dropdown.grid(row=3, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))

tk.Label(frame, text="Contour Intervals (ft or m):").grid(row=4, column=0, sticky=tk.W, pady=(5, 5))
contour_intervals_entry = tk.Entry(frame)
contour_intervals_entry.grid(row=4, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))

tk.Label(frame, text="Scale (ft=0.3048, m=1):").grid(row=5, column=0, sticky=tk.W, pady=(5, 5))
scale_entry = tk.Entry(frame)
scale_entry.grid(row=5, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))

tk.Label(frame, text="Friction Value:").grid(row=6, column=0, sticky=tk.W, pady=(5, 5))
friction_entry = tk.Entry(frame)
friction_entry.grid(row=6, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))
friction_entry.insert(0, '0.044')  # Set default value to 0.04

tk.Label(frame, text="Volume Calculation Method:").grid(row=7, column=0, sticky=tk.W, pady=(5, 5))
volume_method_var = tk.StringVar(root)
volume_method_var.set("Trapezoidal")  # set default value to Trapzoidal

volume_method_dropdown = tk.OptionMenu(frame, volume_method_var, "Trapezoidal", "Prismoidal")
volume_method_dropdown.grid(row=7, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))


tk.Label(frame, text="Waterbody Name:").grid(row=8, column=0, sticky=tk.W, pady=(5, 5))
name_entry = tk.Entry(frame)
name_entry.grid(row=8, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))

tk.Label(frame, text="Output Folder:").grid(row=9, column=0, sticky=tk.W, pady=(5, 5))
output_entry = tk.Entry(frame)
output_entry.grid(row=9, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))
tk.Button(frame, text="Browse", command=lambda: browse_directory(output_entry)).grid(row=9, column=2, sticky=tk.EW, pady=(5, 5))

contour_intervals_entry.insert(0, '5')
scale_entry.insert(0, '0.3048')

# Create a new frame for the second page
frame2 = tk.Frame(root, padx=10, pady=10)

# Function to switch to the second page
# def switch_to_page2():
#     frame.pack_forget()  # Hide the first page
#     frame2.pack(fill=tk.BOTH, expand=True)  # Show the second page

tk.Button(frame, text="Run", command=run_program).grid(row=10, column=1, sticky=tk.EW, padx=(0, 10), pady=(5, 5))
def show_help():
    help_text = """
    This program generates bathymetry for the CE-QUAL-W2 model using Trapezoidal and Prismoidal methods. To run the program, users need to follow these steps:

     1. Provide topographic contours, the river centerline (for each branch), and W2 segments.
     2. Ensure that the topographic lines in the contour shapefile are free of geometry errors such as dangling lines, self-intersecting lines, etc.
     3. Ensure that no W2 segments are intersected by more than one river centerline.
     4. The contour shapefile must have an Elevation attribute, which will be used to generate the bathymetry.
     5. Both the W2 segments and centerline shapefiles must include a 'Segment' column.
     6. Use the W2 active segment numbers for the centerline and W2 segments shapefiles, e.g., the first segment should be 2. The program will add inactive cells to the output automatically.
       Note: In CE-QUAL-W2, each branch or waterbody starts and ends with inactive segments. For example, there are two inactive segments when a branch is added to a waterbody. 
       When adding a branch to the main waterbody, users should account for these inactive segments and name the segments appropriately.

       For inquiries or to report bugs, please contact afshin.shabani64@gamil.com
    """
    help_window = tk.Toplevel(root)
    help_window.title("Notes")
    
    # Add the Tetra-Tech-Logo.png to the help window
    # logo = tk.PhotoImage(file="Tetra-Tech-Logo.png")
    # logo_label = tk.Label(help_window, image=logo)
    # logo_label.image = logo  # Keep a reference to avoid garbage collection
    help_window.iconphoto(False, custom_icon)
    # logo_label.pack(pady=(10, 0))
    help_label = tk.Label(help_window, text=help_text, justify=tk.LEFT, padx=10, pady=10)
    help_label.pack()

# Function to switch to the bathymetry page
def show_bathymetry_page():
    frame.pack(fill=tk.BOTH, expand=True)  # Show the bathymetry page
    frame2.pack_forget()  # Hide the second page
# Create a menu bar
menu_bar = tk.Menu(root)

# Create a Help menu and add it to the menu bar
# Create a Bathymetry menu and add it to the menu bar
# bathymetry_menu = tk.Menu(menu_bar, tearoff=0)
# bathymetry_menu.add_command(label="Bathymetry", command=show_bathymetry_page)
# menu_bar.add_cascade(label="Bathymetry", menu=bathymetry_menu)

# Create a new menu item for the second page
# calibration_menu = tk.Menu(menu_bar, tearoff=0)
# calibration_menu.add_command(label="Bathymetry Calibration" )#command=switch_to_page2
# menu_bar.add_cascade(label="Calibration", menu=calibration_menu)

help_menu = tk.Menu(menu_bar, tearoff=0)
help_menu.add_command(label="Users' Guide", command=show_help)
menu_bar.add_cascade(label="Help", menu=help_menu)

# Display the menu bar
root.config(menu=menu_bar)


root.mainloop()
