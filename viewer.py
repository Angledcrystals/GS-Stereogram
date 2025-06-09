#!/usr/bin/env python3
"""
G-S Divine Stereo Viewer - Complete Implementation
S-Coordinate Pixel Divination System

Author: Angledcrystals
Date: 2025-06-09
Time: 07:15:46 UTC

Complete implementation of G-S coordinate-based stereo visualization
with S-coordinate pixel divination instead of traditional hole filling.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import cv2
from PIL import Image, ImageTk
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import json
import os
import time
from datetime import datetime
from scipy import ndimage
from scipy.spatial import cKDTree
from scipy.interpolate import griddata, RectBivariateSpline

class GSDivineStereoViewer:
    def __init__(self, root):
        self.root = root
        self.root.title("G-S Divine Stereo Viewer v1.3 - S-Coordinate Pixel Divination")
        self.root.geometry("1600x1000")
        
        # Data storage
        self.texture_image = None
        self.depth_map = None
        self.gs_alignment_data = None
        self.gs_coordinate_map = None
        self.subpixel_enhancement_map = None
        self.last_processing_info = None
        self.last_mapping_info = None
        
        # Processed results
        self.left_eye_image = None
        self.right_eye_image = None
        self.stereo_combined = None
        
        # G-S Divine Pixel parameters
        self.gs_divine_enabled = tk.BooleanVar(value=True)
        self.divine_method = tk.StringVar(value="s_coordinate_flow")
        self.divine_precision = tk.DoubleVar(value=0.1)
        self.divine_search_radius = tk.IntVar(value=8)
        self.divine_similarity_threshold = tk.DoubleVar(value=10.0)
        
        # G-S Sub-pixel parameters
        self.gs_subpixel_enabled = tk.BooleanVar(value=True)
        self.gs_enhancement_factor = tk.DoubleVar(value=2.0)
        self.gs_interpolation_method = tk.StringVar(value="gs_geometric")
        self.gs_coordinate_precision = tk.DoubleVar(value=0.1)
        self.use_gs_alignment_data = tk.BooleanVar(value=False)
        
        # Stereo parameters
        self.depth_scale = tk.DoubleVar(value=1.0)
        self.eye_separation = tk.DoubleVar(value=10.0)
        self.depth_offset = tk.DoubleVar(value=0.0)
        self.max_disparity = tk.DoubleVar(value=30.0)
        self.depth_invert = tk.BooleanVar(value=False)
        self.depth_normalize = tk.BooleanVar(value=True)
        
        # Enhanced processing parameters
        self.displacement_method = tk.StringVar(value="gs_divine_planar")
        self.hole_fill_method = tk.StringVar(value="s_coordinate_divine")
        self.edge_handling = tk.StringVar(value="gs_spherical")
        self.anti_aliasing = tk.BooleanVar(value=True)
        self.subpixel_accuracy = tk.BooleanVar(value=True)
        
        # Display parameters
        self.stereo_mode = tk.StringVar(value="side_by_side")
        self.show_guides = tk.BooleanVar(value=True)
        self.show_gs_overlay = tk.BooleanVar(value=False)
        self.show_divine_overlay = tk.BooleanVar(value=False)
        
        # Processing parameters
        self.smoothing_enabled = tk.BooleanVar(value=True)
        self.smoothing_sigma = tk.DoubleVar(value=1.0)
        self.edge_blending = tk.BooleanVar(value=True)
        self.quality_mode = tk.StringVar(value="gs_divine_enhanced")
        
        # Debug and visualization
        self.debug_divine_process = tk.BooleanVar(value=False)
        self.show_divine_statistics = tk.BooleanVar(value=True)
        self.save_divine_maps = tk.BooleanVar(value=False)
        
        self.setup_gui()
        
    def setup_gui(self):
        """Setup the main GUI layout."""
        control_frame = ttk.Frame(self.root)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        viz_frame = ttk.Frame(self.root)
        viz_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        self.setup_scrollable_controls(control_frame)
        self.setup_stereo_visualization(viz_frame)
        
    def setup_scrollable_controls(self, parent):
        """Setup scrollable control panel."""
        canvas = tk.Canvas(parent, highlightthickness=0, width=480)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Mouse wheel scrolling
        def _on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        def _bind_mousewheel(event):
            canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        def _unbind_mousewheel(event):
            canvas.unbind_all("<MouseWheel>")
        
        scrollable_frame.bind('<Enter>', _bind_mousewheel)
        scrollable_frame.bind('<Leave>', _unbind_mousewheel)
        canvas.bind('<Enter>', _bind_mousewheel)
        canvas.bind('<Leave>', _unbind_mousewheel)
        
        self.setup_control_panel(scrollable_frame)
        
    def setup_control_panel(self, parent):
        """Setup the control panel with G-S divine pixel options."""
        
        # Title
        title_label = ttk.Label(parent, text="G-S Divine Stereo Viewer", font=("Arial", 16, "bold"))
        title_label.pack(pady=(0, 5))
        
        subtitle_label = ttk.Label(parent, text="S-Coordinate Pixel Divination v1.3", font=("Arial", 10))
        subtitle_label.pack(pady=(0, 15))
        
        # Data Loading Section
        data_frame = ttk.LabelFrame(parent, text="Step 1: Load G-S Data", padding=10)
        data_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Button(data_frame, text="🖼️ Load Texture Image", 
                  command=self.load_texture_image, width=40).pack(fill=tk.X, pady=2)
        
        ttk.Button(data_frame, text="🗺️ Load G-S Depth Map", 
                  command=self.load_depth_map, width=40).pack(fill=tk.X, pady=2)
        
        ttk.Button(data_frame, text="📊 Load G-S Alignment Data (Optional)", 
                  command=self.load_gs_alignment_data, width=40).pack(fill=tk.X, pady=2)
        
        ttk.Button(data_frame, text="🔄 Load All (Auto-detect)", 
                  command=self.load_all_auto, width=40).pack(fill=tk.X, pady=2)
        
        self.data_status_label = ttk.Label(data_frame, text="No data loaded", 
                                          foreground="red")
        self.data_status_label.pack(pady=5)
        
        # G-S Divine Pixel Section (NEW MAIN FEATURE)
        divine_frame = ttk.LabelFrame(parent, text="Step 2: S-Coordinate Pixel Divination", padding=10)
        divine_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(divine_frame, text="Enable S-coordinate pixel divination", 
                       variable=self.gs_divine_enabled, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        ttk.Label(divine_frame, text="Divine Method:").pack(anchor=tk.W)
        divine_method_combo = ttk.Combobox(divine_frame, textvariable=self.divine_method,
                                          values=["s_coordinate_flow", "s_coordinate_similarity", 
                                                 "geometric_extrapolation", "combined_divine", "gs_enhanced_divine"])
        divine_method_combo.pack(fill=tk.X, pady=2)
        divine_method_combo.bind('<<ComboboxSelected>>', self.on_parameter_change)
        
        ttk.Label(divine_frame, text="Divine Precision (S-coordinate sensitivity):").pack(anchor=tk.W)
        divine_precision_frame = ttk.Frame(divine_frame)
        divine_precision_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(divine_precision_frame, from_=0.01, to=1.0, variable=self.divine_precision, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(divine_precision_frame, textvariable=self.divine_precision, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(divine_frame, text="Search Radius (pixels):").pack(anchor=tk.W)
        search_radius_frame = ttk.Frame(divine_frame)
        search_radius_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(search_radius_frame, from_=3, to=20, variable=self.divine_search_radius, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(search_radius_frame, textvariable=self.divine_search_radius, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(divine_frame, text="Similarity Threshold (degrees):").pack(anchor=tk.W)
        similarity_frame = ttk.Frame(divine_frame)
        similarity_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(similarity_frame, from_=1.0, to=30.0, variable=self.divine_similarity_threshold, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(similarity_frame, textvariable=self.divine_similarity_threshold, width=8).pack(side=tk.RIGHT)
        
        # Debug options for divine process
        ttk.Checkbutton(divine_frame, text="Debug divine process", 
                       variable=self.debug_divine_process, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        ttk.Checkbutton(divine_frame, text="Show divine statistics", 
                       variable=self.show_divine_statistics, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        self.test_divine_button = ttk.Button(divine_frame, text="🔮 Test Divine Process", 
                                           command=self.test_divine_process, 
                                           state="disabled")
        self.test_divine_button.pack(fill=tk.X, pady=5)
        
        # G-S Sub-pixel Enhancement Section  
        gs_frame = ttk.LabelFrame(parent, text="Step 3: G-S Sub-pixel Enhancement", padding=10)
        gs_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(gs_frame, text="Enable G-S sub-pixel enhancement", 
                       variable=self.gs_subpixel_enabled, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        ttk.Label(gs_frame, text="Enhancement Factor (sub-pixel resolution):").pack(anchor=tk.W)
        enhance_frame = ttk.Frame(gs_frame)
        enhance_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(enhance_frame, from_=1.0, to=8.0, variable=self.gs_enhancement_factor, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(enhance_frame, textvariable=self.gs_enhancement_factor, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(gs_frame, text="G-S Interpolation Method:").pack(anchor=tk.W)
        gs_interp_combo = ttk.Combobox(gs_frame, textvariable=self.gs_interpolation_method,
                                      values=["gs_geometric", "gs_spherical", "gs_coordinate", 
                                             "bilinear", "bicubic", "gs_hadit_enhanced"])
        gs_interp_combo.pack(fill=tk.X, pady=2)
        gs_interp_combo.bind('<<ComboboxSelected>>', self.on_parameter_change)
        
        ttk.Label(gs_frame, text="G-S Coordinate Precision:").pack(anchor=tk.W)
        precision_frame = ttk.Frame(gs_frame)
        precision_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(precision_frame, from_=0.01, to=1.0, variable=self.gs_coordinate_precision, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(precision_frame, textvariable=self.gs_coordinate_precision, width=8).pack(side=tk.RIGHT)
        
        ttk.Checkbutton(gs_frame, text="Use G-S alignment data for enhancement", 
                       variable=self.use_gs_alignment_data, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        self.generate_gs_coords_button = ttk.Button(gs_frame, text="🧮 Generate G-S Coordinate Maps", 
                                                   command=self.generate_gs_coordinate_maps, 
                                                   state="disabled")
        self.generate_gs_coords_button.pack(fill=tk.X, pady=5)
        
        # Depth Processing Section
        depth_frame = ttk.LabelFrame(parent, text="Step 4: Depth Map Processing", padding=10)
        depth_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(depth_frame, text="Depth Scale Factor:").pack(anchor=tk.W)
        depth_scale_frame = ttk.Frame(depth_frame)
        depth_scale_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(depth_scale_frame, from_=0.1, to=3.0, variable=self.depth_scale, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(depth_scale_frame, textvariable=self.depth_scale, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(depth_frame, text="Depth Offset:").pack(anchor=tk.W)
        offset_frame = ttk.Frame(depth_frame)
        offset_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(offset_frame, from_=-50.0, to=50.0, variable=self.depth_offset, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(offset_frame, textvariable=self.depth_offset, width=8).pack(side=tk.RIGHT)
        
        ttk.Checkbutton(depth_frame, text="Invert depth values", 
                       variable=self.depth_invert, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        ttk.Checkbutton(depth_frame, text="Normalize depth to [0,1]", 
                       variable=self.depth_normalize, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        # Enhanced Displacement Section
        enhanced_frame = ttk.LabelFrame(parent, text="Step 5: Divine Enhanced Displacement", padding=10)
        enhanced_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(enhanced_frame, text="Displacement Method:").pack(anchor=tk.W)
        disp_combo = ttk.Combobox(enhanced_frame, textvariable=self.displacement_method,
                                 values=["gs_divine_planar", "gs_divine_spherical", "gs_enhanced", 
                                        "gs_coordinate_based", "backward_mapping", "forward_mapping"])
        disp_combo.pack(fill=tk.X, pady=2)
        disp_combo.bind('<<ComboboxSelected>>', self.on_parameter_change)
        
        ttk.Label(enhanced_frame, text="Hole Replacement Method:").pack(anchor=tk.W)
        hole_combo = ttk.Combobox(enhanced_frame, textvariable=self.hole_fill_method,
                                 values=["s_coordinate_divine", "gs_flow_divine", "gs_similarity_divine",
                                        "geometric_divine", "inpainting", "interpolation"])
        hole_combo.pack(fill=tk.X, pady=2)
        hole_combo.bind('<<ComboboxSelected>>', self.on_parameter_change)
        
        ttk.Label(enhanced_frame, text="Edge Handling:").pack(anchor=tk.W)
        edge_combo = ttk.Combobox(enhanced_frame, textvariable=self.edge_handling,
                                 values=["gs_spherical", "gs_coordinate", "clamp", "wrap", "mirror"])
        edge_combo.pack(fill=tk.X, pady=2)
        edge_combo.bind('<<ComboboxSelected>>', self.on_parameter_change)
        
        # Stereo Parameters Section
        stereo_frame = ttk.LabelFrame(parent, text="Step 6: Stereo Parameters", padding=10)
        stereo_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Label(stereo_frame, text="Eye Separation (pixels):").pack(anchor=tk.W)
        eye_sep_frame = ttk.Frame(stereo_frame)
        eye_sep_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(eye_sep_frame, from_=1.0, to=30.0, variable=self.eye_separation, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(eye_sep_frame, textvariable=self.eye_separation, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(stereo_frame, text="Maximum Disparity:").pack(anchor=tk.W)
        max_disp_frame = ttk.Frame(stereo_frame)
        max_disp_frame.pack(fill=tk.X, pady=2)
        ttk.Scale(max_disp_frame, from_=5.0, to=100.0, variable=self.max_disparity, 
                 orient=tk.HORIZONTAL, command=self.on_parameter_change).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Label(max_disp_frame, textvariable=self.max_disparity, width=8).pack(side=tk.RIGHT)
        
        ttk.Label(stereo_frame, text="Stereo Mode:").pack(anchor=tk.W, pady=(10, 0))
        mode_combo = ttk.Combobox(stereo_frame, textvariable=self.stereo_mode,
                                 values=["side_by_side", "over_under", "anaglyph_red_cyan", 
                                        "left_only", "right_only", "cross_eyed"])
        mode_combo.pack(fill=tk.X, pady=2)
        mode_combo.bind('<<ComboboxSelected>>', self.on_parameter_change)
        
        # Visualization Options
        viz_options_frame = ttk.LabelFrame(parent, text="Visualization Options", padding=10)
        viz_options_frame.pack(fill=tk.X, pady=(0, 10))
        
        ttk.Checkbutton(viz_options_frame, text="Show alignment guides", 
                       variable=self.show_guides, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        ttk.Checkbutton(viz_options_frame, text="Show G-S coordinate overlay", 
                       variable=self.show_gs_overlay, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        ttk.Checkbutton(viz_options_frame, text="Show divine process overlay", 
                       variable=self.show_divine_overlay, command=self.on_parameter_change).pack(anchor=tk.W, pady=2)
        
        # Processing Section
        process_frame = ttk.LabelFrame(parent, text="Step 7: Generate Divine Stereo", padding=10)
        process_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.generate_button = ttk.Button(process_frame, text="🔮 Generate Divine G-S Stereo", 
                                         command=self.generate_divine_gs_stereo, 
                                         state="disabled")
        self.generate_button.pack(fill=tk.X, pady=5)
        
        self.auto_update_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(process_frame, text="Auto-update on parameter change", 
                       variable=self.auto_update_var).pack(anchor=tk.W, pady=2)
        
        self.process_status_label = ttk.Label(process_frame, text="Load G-S data first", 
                                             foreground="orange")
        self.process_status_label.pack(pady=5)
        
        # Progress bar
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(process_frame, variable=self.progress_var, 
                                          maximum=100)
        self.progress_bar.pack(fill=tk.X, pady=5)
        
        # Export Section
        export_frame = ttk.LabelFrame(parent, text="Step 8: Export Divine Results", padding=10)
        export_frame.pack(fill=tk.X, pady=(0, 10))
        
        self.export_button = ttk.Button(export_frame, text="💾 Export Divine Stereo", 
                                       command=self.export_divine_stereo, 
                                       state="disabled")
        self.export_button.pack(fill=tk.X, pady=2)
        
        self.export_divine_maps_button = ttk.Button(export_frame, text="🗺️ Export Divine Maps", 
                                                   command=self.export_divine_maps, 
                                                   state="disabled")
        self.export_divine_maps_button.pack(fill=tk.X, pady=2)
        
        # Statistics Section
        stats_frame = ttk.LabelFrame(parent, text="Divine Enhancement Statistics", padding=5)
        stats_frame.pack(fill=tk.X, pady=(10, 0))
        
        self.stats_label = ttk.Label(stats_frame, text="No divine enhancement generated", 
                                    font=("Courier", 8))
        self.stats_label.pack(fill=tk.X)
        
    def setup_stereo_visualization(self, parent):
        """Setup the stereo visualization panel."""
        self.fig = Figure(figsize=(18, 14), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, parent)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Toolbar
        toolbar_frame = ttk.Frame(parent)
        toolbar_frame.pack(fill=tk.X)
        
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()
        
        self.show_welcome_display()
        
    def show_welcome_display(self):
        """Show welcome message."""
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        
        welcome_text = """G-S DIVINE STEREO VIEWER v1.3 - S-COORDINATE PIXEL DIVINATION
by Angledcrystals

🔮 Divine Pixel Features:
• S-coordinate flow-based pixel divination
• Geometric extrapolation using G-S relationships  
• Similarity matching in S-coordinate space
• True pixel "divination" instead of hole filling
• Multi-method divine process with fallbacks

🧮 How Divine Process Works:
Instead of filling holes with nearby pixels, the system
DIVINES what pixel should naturally exist at each location
based on your G-S coordinate system's geometric relationships.

🌀 S-Coordinate Flow Tracing:
Traces along S-coordinate flow lines to find the natural
source pixel that should appear at displaced locations.

📐 Geometric Synthesis:
Uses sphere intersection mathematics to synthesize pixels
that should exist based on G-S coordinate relationships.

Load your G-S data to unlock divine pixel enhancement!"""
        
        ax.text(0.5, 0.5, welcome_text, ha='center', va='center', 
                fontsize=10, transform=ax.transAxes, 
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcyan", alpha=0.9))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        ax.set_title("G-S Divine Pixel Enhanced Stereo Viewer", fontsize=16, weight='bold')
        
        self.canvas.draw()
        
    # === DATA LOADING METHODS ===
    
    def load_texture_image(self):
        """Load texture image file."""
        file_path = filedialog.askopenfilename(
            title="Select Texture Image",
            filetypes=[("Image files", "*.png *.jpg *.jpeg *.bmp *.tiff"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                # Load image
                pil_image = Image.open(file_path)
                self.texture_image = np.array(pil_image.convert('RGB'))
                
                print(f"Loaded texture image: {self.texture_image.shape}")
                self.update_data_status()
                self.show_loaded_data()
                
                # Enable coordinate generation if depth is also loaded
                if self.depth_map is not None:
                    self.generate_gs_coords_button.config(state="normal")
                    
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load texture image: {str(e)}")
                
    def load_depth_map(self):
        """Load G-S depth map file."""
        file_path = filedialog.askopenfilename(
            title="Select G-S Depth Map",
            filetypes=[("Image files", "*.png *.jpg *.jpeg *.bmp *.tiff"), 
                      ("NumPy files", "*.npy"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                if file_path.endswith('.npy'):
                    self.depth_map = np.load(file_path)
                else:
                    pil_image = Image.open(file_path)
                    if pil_image.mode == 'RGB':
                        # Convert RGB to grayscale
                        gray_image = pil_image.convert('L')
                        self.depth_map = np.array(gray_image) / 255.0
                    else:
                        self.depth_map = np.array(pil_image) / 255.0
                
                print(f"Loaded G-S depth map: {self.depth_map.shape}")
                print(f"Depth range: [{self.depth_map.min():.3f}, {self.depth_map.max():.3f}]")
                
                self.update_data_status()
                self.show_loaded_data()
                
                # Enable coordinate generation if texture is also loaded
                if self.texture_image is not None:
                    self.generate_gs_coords_button.config(state="normal")
                    
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load depth map: {str(e)}")

    def load_gs_alignment_data(self):
        """Load S-NUIT alignment data from JSON file."""
        file_path = filedialog.askopenfilename(
            title="Load S-NUIT Alignment JSON",
            filetypes=[
                ("JSON files", "*.json"),
                ("All files", "*.*")
            ]
        )

        if not file_path:
            return

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                alignment_json = json.load(f)

            # Extract alignment data
            alignments = alignment_json.get('alignments', [])
            metadata = alignment_json.get('metadata', {})
            condition_summary = alignment_json.get('condition_summary', {})
            statistics = alignment_json.get('statistics', {})

            # Check if data is organized by G-S coordinates
            if metadata.get('organized_by_gs_coordinates', False):
                # Flatten organized data
                flattened_alignments = []
                for condition, group_alignments in alignments.items():
                    flattened_alignments.extend(group_alignments)
                alignments = flattened_alignments

            # Validate required fields
            if alignments:
                required_fields = ['G_theta_deg', 'G_phi_deg', 'S_x', 'S_y']
                sample = alignments[0]
                missing_fields = [field for field in required_fields if field not in sample]
                if missing_fields:
                    raise ValueError(f"Alignment data missing required fields: {missing_fields}")

            # Store alignment data
            self.gs_alignment_data = alignments
            self.gs_alignment_metadata = {
                'source': metadata.get('source_file', 'Unknown'),
                'generated': metadata.get('generated', 'Unknown'),
                'total_alignments': len(alignments),
                'total_conditions': metadata.get('total_conditions', 0),
                'condition_summary': condition_summary,
                'statistics': statistics
            }

            # Update UI
            alignment_count = len(alignments)
            condition_count = len(condition_summary)

            self.process_status_label.config(
                text=f"✅ Loaded {alignment_count:,} S-NUIT alignments ({condition_count} conditions)",
                foreground="green"
            )

            # Enable coordinate map generation
            if self.texture_image is not None and self.depth_map is not None:
                self.generate_gs_coords_button.config(state="normal")

            # Update data status
            self.update_data_status()

            # Update statistics
            self.update_statistics()

            print(f"🔮 Loaded S-NUIT alignment data:")
            print(f"   Alignments: {alignment_count:,}")
            print(f"   Conditions: {condition_count}")
            print(f"   Source: {metadata.get('source_file', 'Unknown')}")
            print(f"   Generated: {metadata.get('generated', 'Unknown')}")

            if statistics:
                print(f"   G_theta range: {statistics.get('g_theta_range', 'N/A')}")
                print(f"   G_phi range: {statistics.get('g_phi_range', 'N/A')}")
                print(f"   S_x range: {statistics.get('s_x_range', 'N/A')}")
                print(f"   S_y range: {statistics.get('s_y_range', 'N/A')}")

            messagebox.showinfo(
                "Alignment Data Loaded",
                f"Successfully loaded {alignment_count:,} S-NUIT alignments\n"
                f"Conditions analyzed: {condition_count}\n"
                f"Source: {metadata.get('source_file', 'Unknown')}\n"
                f"Generated: {metadata.get('generated', 'Unknown')}"
            )

        except Exception as e:
            messagebox.showerror("Load Error", f"Failed to load alignment JSON: {str(e)}")
            print(f"Error loading alignment JSON: {e}")
                
    def load_all_auto(self):
        """Auto-detect and load G-S data files from directory."""
        directory = filedialog.askdirectory(title="Select directory with G-S data")
        
        if directory:
            try:
                loaded_files = []
                
                # Look for common file patterns
                for filename in os.listdir(directory):
                    filepath = os.path.join(directory, filename)
                    
                    # Texture patterns
                    if any(pattern in filename.lower() for pattern in ['texture', 'color', 'rgb', 'image']):
                        if filename.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tiff')):
                            pil_image = Image.open(filepath)
                            self.texture_image = np.array(pil_image.convert('RGB'))
                            loaded_files.append(f"Texture: {filename}")
                    
                    # Depth map patterns
                    elif any(pattern in filename.lower() for pattern in ['depth', 'height', 'displacement']):
                        if filename.lower().endswith('.npy'):
                            self.depth_map = np.load(filepath)
                            loaded_files.append(f"Depth: {filename}")
                        elif filename.lower().endswith(('.png', '.jpg', '.jpeg', '.bmp', '.tiff')):
                            pil_image = Image.open(filepath)
                            if pil_image.mode == 'RGB':
                                gray_image = pil_image.convert('L')
                                self.depth_map = np.array(gray_image) / 255.0
                            else:
                                self.depth_map = np.array(pil_image) / 255.0
                            loaded_files.append(f"Depth: {filename}")
                    
                    # Alignment data patterns
                    elif any(pattern in filename.lower() for pattern in ['alignment', 'gs', 'coord']):
                        if filename.lower().endswith('.json'):
                            with open(filepath, 'r') as f:
                                self.gs_alignment_data = json.load(f)
                            loaded_files.append(f"Alignments: {filename}")
                
                if loaded_files:
                    self.update_data_status()
                    self.show_loaded_data()
                    
                    if self.texture_image is not None and self.depth_map is not None:
                        self.generate_gs_coords_button.config(state="normal")
                    
                    messagebox.showinfo("Auto-load Success", 
                                      f"Loaded files:\n" + "\n".join(loaded_files))
                else:
                    messagebox.showwarning("Auto-load", "No G-S data files found in directory")
                    
            except Exception as e:
                messagebox.showerror("Auto-load Error", f"Failed to auto-load: {str(e)}")
                
    # === DEPTH PROCESSING METHODS ===
    
    def apply_depth_processing(self, depth_map):
        """Apply depth processing parameters."""
        if self.depth_invert.get():
            depth_map = 1.0 - depth_map
        
        depth_map = depth_map * self.depth_scale.get()
        depth_map = depth_map + self.depth_offset.get() / 100.0
        
        if self.smoothing_enabled.get() and self.smoothing_sigma.get() > 0:
            depth_map = ndimage.gaussian_filter(depth_map, sigma=self.smoothing_sigma.get())
        
        if self.depth_normalize.get():
            if depth_map.max() > depth_map.min():
                depth_map = (depth_map - depth_map.min()) / (depth_map.max() - depth_map.min())
        
        depth_map = np.clip(depth_map, 0.0, 1.0)
        return depth_map
        
    def generate_gs_coordinate_maps(self):
        """Generate G-S coordinate maps from depth map and alignment data."""
        if self.depth_map is None:
            messagebox.showwarning("No Data", "Load G-S depth map first")
            return
            
        try:
            self.progress_var.set(0)
            start_time = time.time()
            
            height, width = self.depth_map.shape
            
            # Initialize coordinate maps
            g_theta_map = np.zeros((height, width))
            g_phi_map = np.zeros((height, width))
            s_x_map = np.zeros((height, width))
            s_y_map = np.zeros((height, width))
            
            self.progress_var.set(20)
            self.root.update_idletasks()
            
            # Method 1: If we have alignment data, use it for ground truth
            if self.gs_alignment_data and self.use_gs_alignment_data.get():
                g_theta_map, g_phi_map, s_x_map, s_y_map = self.interpolate_from_alignment_data(height, width)
            else:
                # Method 2: Derive G-S coordinates from depth map geometry
                g_theta_map, g_phi_map, s_x_map, s_y_map = self.derive_gs_coordinates_from_depth(height, width)
            
            self.progress_var.set(60)
            self.root.update_idletasks()
            
            # Store coordinate maps
            self.gs_coordinate_map = {
                'g_theta': g_theta_map,
                'g_phi': g_phi_map,
                's_x': s_x_map,
                's_y': s_y_map
            }
            
            # Generate sub-pixel enhancement map
            if self.gs_subpixel_enabled.get():
                self.subpixel_enhancement_map = self.generate_subpixel_enhancement(
                    g_theta_map, g_phi_map, s_x_map, s_y_map
                )
            
            self.progress_var.set(80)
            self.root.update_idletasks()
            
            processing_time = time.time() - start_time
            
            self.process_status_label.config(
                text=f"✅ G-S coordinate maps generated ({processing_time:.2f}s)", 
                foreground="green"
            )
            
            # Enable divine process testing
            self.test_divine_button.config(state="normal")
            
            print(f"Generated G-S coordinate maps in {processing_time:.2f} seconds")
            self.visualize_gs_coordinate_maps()
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate G-S coordinate maps: {str(e)}")
        finally:
            self.progress_var.set(100)

    def load_gs_alignment_data(self):
        """Load G-S alignment data file with proper JSON structure handling."""
        file_path = filedialog.askopenfilename(
            title="Select G-S Alignment Data",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )

        if file_path:
            try:
                with open(file_path, 'r') as f:
                    alignment_json = json.load(f)

                # Handle the converted JSON structure
                if 'alignments' in alignment_json:
                    # This is the new JSON format from the converter
                    alignments = alignment_json['alignments']
                    metadata = alignment_json.get('metadata', {})

                    # Check if data is organized by G-S coordinates
                    if metadata.get('organized_by_gs_coordinates', False):
                        # Flatten organized data
                        flattened_alignments = []
                        for condition, group_alignments in alignments.items():
                            flattened_alignments.extend(group_alignments)
                        self.gs_alignment_data = flattened_alignments
                    else:
                        self.gs_alignment_data = alignments

                    # Store metadata
                    self.gs_alignment_metadata = metadata

                    print(f"Loaded G-S alignment data: {len(self.gs_alignment_data)} alignments")
                    print(f"Metadata: {metadata.get('source_file', 'Unknown source')}")

                    # Print sample alignment for verification
                    if self.gs_alignment_data:
                        sample = self.gs_alignment_data[0]
                        print(f"Sample alignment keys: {list(sample.keys())}")

                        # Verify required fields exist
                        required_fields = ['G_theta_deg', 'G_phi_deg', 'S_x', 'S_y']
                        missing_fields = [field for field in required_fields if field not in sample]
                        if missing_fields:
                            print(f"Warning: Missing required fields: {missing_fields}")

                else:
                    # This is the old format (direct list of alignments)
                    self.gs_alignment_data = alignment_json
                    print(f"Loaded G-S alignment data (old format): {len(self.gs_alignment_data)} alignments")

                self.update_data_status()

            except Exception as e:
                messagebox.showerror("Error", f"Failed to load alignment data: {str(e)}")
                print(f"Error details: {str(e)}")

    def derive_gs_coordinates_from_depth(self, height, width):
        """Derive G-S coordinates from depth map geometry."""
        # Create coordinate grids
        y_grid, x_grid = np.mgrid[0:height, 0:width]
        
        # Normalize coordinates to [-1, 1] range
        x_normalized = (x_grid / width) * 2.0 - 1.0
        y_normalized = (y_grid / height) * 2.0 - 1.0
        
        # Use depth to estimate spherical coordinates
        depth_values = self.depth_map
        
        # Estimate G_theta and G_phi from image coordinates and depth
        g_theta_est = np.arctan2(y_normalized, x_normalized) * 180.0 / np.pi
        
        # Estimate phi based on distance from center and depth
        radius_from_center = np.sqrt(x_normalized**2 + y_normalized**2)
        g_phi_est = (radius_from_center * depth_values) * 180.0 / np.pi
        
        # Clamp phi to valid range
        g_phi_est = np.clip(g_phi_est, 0, 180)
        
        # S coordinates are approximated as normalized image coordinates
        s_x_est = x_normalized
        s_y_est = y_normalized
        
        return g_theta_est, g_phi_est, s_x_est, s_y_est
        
    def generate_subpixel_enhancement(self, g_theta, g_phi, s_x, s_y):
        """Generate sub-pixel enhancement map using G-S coordinate interpolation."""
        enhancement_factor = int(self.gs_enhancement_factor.get())
        if enhancement_factor <= 1:
            return None
        
        height, width = g_theta.shape
        new_height = height * enhancement_factor
        new_width = width * enhancement_factor
        
        # Create high-resolution coordinate grids
        y_hr, x_hr = np.mgrid[0:new_height, 0:new_width]
        
        # Map high-res coordinates back to original coordinate space
        y_orig = y_hr / enhancement_factor
        x_orig = x_hr / enhancement_factor
        
        method = self.gs_interpolation_method.get()
        
        if method == "gs_geometric":
            enhanced_map = self.gs_geometric_interpolation(g_theta, g_phi, s_x, s_y, x_orig, y_orig)
        elif method == "gs_spherical":
            enhanced_map = self.gs_spherical_interpolation(g_theta, g_phi, x_orig, y_orig)
        elif method == "gs_coordinate":
            enhanced_map = self.gs_coordinate_interpolation(s_x, s_y, x_orig, y_orig)
        elif method == "gs_hadit_enhanced":
            enhanced_map = self.gs_hadit_interpolation(g_theta, g_phi, s_x, s_y, x_orig, y_orig)
        else:
            enhanced_map = self.standard_interpolation(self.depth_map, x_orig, y_orig, method)
        
        return enhanced_map
        
    def gs_geometric_interpolation(self, g_theta, g_phi, s_x, s_y, x_hr, y_hr):
        """Interpolate using G-S geometric relationships."""
        # Convert G-S spherical coordinates to 3D vectors
        g_theta_rad = g_theta * np.pi / 180.0
        g_phi_rad = g_phi * np.pi / 180.0
        
        # Calculate 3D G vectors
        g_x = np.sin(g_phi_rad) * np.cos(g_theta_rad)
        g_y = np.sin(g_phi_rad) * np.sin(g_theta_rad)
        g_z = np.cos(g_phi_rad)
        
        # Interpolate each component using spline interpolation
        height, width = g_theta.shape
        
        # Create spline interpolators for each component
        y_coords = np.arange(height)
        x_coords = np.arange(width)
        
        spline_gx = RectBivariateSpline(y_coords, x_coords, g_x, kx=3, ky=3)
        spline_gy = RectBivariateSpline(y_coords, x_coords, g_y, kx=3, ky=3)
        spline_gz = RectBivariateSpline(y_coords, x_coords, g_z, kx=3, ky=3)
        
        # Evaluate at high-resolution coordinates
        g_x_hr = spline_gx.ev(y_hr, x_hr)
        g_y_hr = spline_gy.ev(y_hr, x_hr)
        g_z_hr = spline_gz.ev(y_hr, x_hr)
        
        # Convert back to depth-like values
        enhanced_depth = np.sqrt(g_x_hr**2 + g_y_hr**2 + g_z_hr**2)
        
        return enhanced_depth
        
    def gs_spherical_interpolation(self, g_theta, g_phi, x_hr, y_hr):
        """Interpolate using spherical coordinate smoothness."""
        height, width = g_theta.shape
        
        # Create spline interpolators for spherical coordinates
        y_coords = np.arange(height)
        x_coords = np.arange(width)
        
        # Handle theta wraparound (angular coordinate)
        g_theta_cos = np.cos(g_theta * np.pi / 180.0)
        g_theta_sin = np.sin(g_theta * np.pi / 180.0)
        
        spline_theta_cos = RectBivariateSpline(y_coords, x_coords, g_theta_cos, kx=3, ky=3)
        spline_theta_sin = RectBivariateSpline(y_coords, x_coords, g_theta_sin, kx=3, ky=3)
        spline_phi = RectBivariateSpline(y_coords, x_coords, g_phi, kx=3, ky=3)
        
        # Evaluate at high-resolution
        theta_cos_hr = spline_theta_cos.ev(y_hr, x_hr)
        theta_sin_hr = spline_theta_sin.ev(y_hr, x_hr)
        phi_hr = spline_phi.ev(y_hr, x_hr)
        
        # Reconstruct theta
        theta_hr = np.arctan2(theta_sin_hr, theta_cos_hr) * 180.0 / np.pi
        
        # Convert to depth using spherical relationship
        enhanced_depth = np.sin(phi_hr * np.pi / 180.0)
        
        return enhanced_depth
        
    def gs_coordinate_interpolation(self, s_x, s_y, x_hr, y_hr):
        """Interpolate using S-coordinate relationships."""
        height, width = s_x.shape
        
        # Create spline interpolators for S coordinates
        y_coords = np.arange(height)
        x_coords = np.arange(width)
        
        spline_sx = RectBivariateSpline(y_coords, x_coords, s_x, kx=3, ky=3)
        spline_sy = RectBivariateSpline(y_coords, x_coords, s_y, kx=3, ky=3)
        
        # Evaluate at high-resolution
        s_x_hr = spline_sx.ev(y_hr, x_hr)
        s_y_hr = spline_sy.ev(y_hr, x_hr)
        
        # Calculate depth from S-coordinate geometry
        radius_sq = 1.0  # Normalized sphere
        s_magnitude_sq = s_x_hr**2 + s_y_hr**2
        
        # Clamp to valid range
        s_magnitude_sq = np.clip(s_magnitude_sq, 0, radius_sq)
        enhanced_depth = np.sqrt(radius_sq - s_magnitude_sq)
        
        return enhanced_depth

    def gs_hadit_interpolation(self, g_theta, g_phi, s_x, s_y, x_hr, y_hr):
        """Enhanced interpolation using Hadit vector relationships."""
        if not self.gs_alignment_data:
            return self.gs_geometric_interpolation(g_theta, g_phi, s_x, s_y, x_hr, y_hr)

        try:
            # Extract Hadit information from alignment data
            hadit_thetas = []
            hadit_phis = []

            for alignment in self.gs_alignment_data:
                # Handle both possible field name formats
                hadit_theta = alignment.get('Hadit_theta_deg',
                                            alignment.get('Hadit_theta', alignment.get('hadit_theta', None)))
                hadit_phi = alignment.get('Hadit_phi_deg', alignment.get('Hadit_phi', alignment.get('hadit_phi', None)))

                if hadit_theta is not None and hadit_phi is not None:
                    hadit_thetas.append(float(hadit_theta))
                    hadit_phis.append(float(hadit_phi))

            if len(hadit_thetas) == 0:
                print("No Hadit data found, falling back to geometric interpolation")
                return self.gs_geometric_interpolation(g_theta, g_phi, s_x, s_y, x_hr, y_hr)

            # Use average Hadit vector for enhancement
            avg_hadit_theta = np.mean(hadit_thetas) * np.pi / 180.0
            avg_hadit_phi = np.mean(hadit_phis) * np.pi / 180.0

            print(f"Using Hadit enhancement with {len(hadit_thetas)} Hadit vectors")
            print(f"Average Hadit: theta={np.mean(hadit_thetas):.1f}°, phi={np.mean(hadit_phis):.1f}°")

            # Calculate Hadit 3D vector
            hadit_x = np.sin(avg_hadit_phi) * np.cos(avg_hadit_theta)
            hadit_y = np.sin(avg_hadit_phi) * np.sin(avg_hadit_theta)
            hadit_z = np.cos(avg_hadit_phi)

            # Use Hadit vector to modify interpolation
            enhanced_depth = self.gs_geometric_interpolation(g_theta, g_phi, s_x, s_y, x_hr, y_hr)

            # Apply Hadit-based modulation
            y_hr_norm = (y_hr - y_hr.mean()) / (y_hr.std() + 1e-6)
            x_hr_norm = (x_hr - x_hr.mean()) / (x_hr.std() + 1e-6)

            hadit_modulation = (hadit_x * x_hr_norm + hadit_y * y_hr_norm) * 0.1
            enhanced_depth = enhanced_depth + hadit_modulation

            return np.clip(enhanced_depth, 0, 1)

        except Exception as e:
            print(f"Error in Hadit interpolation: {e}")
            return self.gs_geometric_interpolation(g_theta, g_phi, s_x, s_y, x_hr, y_hr)
        
    def standard_interpolation(self, depth_map, x_hr, y_hr, method):
        """Standard interpolation methods as fallback."""
        height, width = depth_map.shape
        y_coords = np.arange(height)
        x_coords = np.arange(width)
        
        if method == "bicubic":
            spline = RectBivariateSpline(y_coords, x_coords, depth_map, kx=3, ky=3)
        else:  # bilinear
            spline = RectBivariateSpline(y_coords, x_coords, depth_map, kx=1, ky=1)
        
        return spline.ev(y_hr, x_hr)
            
    # === S-COORDINATE DIVINE PIXEL SYSTEM ===
    
    def divine_pixels_from_s_coordinates(self, displaced_image, original_texture, s_x_map, s_y_map, holes_mask):
        """
        Divine pixels using S-coordinate geometric relationships instead of hole filling.
        
        Args:
            displaced_image: Image after displacement with holes
            original_texture: Original texture image
            s_x_map: S-coordinate X values
            s_y_map: S-coordinate Y values  
            holes_mask: Boolean mask of hole locations
        
        Returns:
            divined_image: Image with geometrically divined pixels
        """
        if self.debug_divine_process.get():
            print("🔮 Divining pixels using S-coordinate geometry...")
        
        height, width = displaced_image.shape[:2]
        divined_image = displaced_image.copy()
        
        # Count initial holes
        initial_holes = np.sum(holes_mask)
        if self.debug_divine_process.get():
            print(f"  Initial holes to divine: {initial_holes}")
        
        # Method selection based on user choice
        divine_method = self.divine_method.get()
        
        if divine_method == "s_coordinate_flow":
            divined_image = self.s_coordinate_flow_divination(
                divined_image, original_texture, s_x_map, s_y_map, holes_mask)
        elif divine_method == "s_coordinate_similarity":
            divined_image = self.s_coordinate_similarity_divination(
                divined_image, original_texture, s_x_map, s_y_map, holes_mask)
        elif divine_method == "geometric_extrapolation":
            divined_image = self.s_coordinate_geometric_extrapolation(
                divined_image, original_texture, s_x_map, s_y_map, holes_mask)
        elif divine_method == "combined_divine":
            # Multi-method approach: try each method in sequence
            divined_image = self.combined_divine_process(
                divined_image, original_texture, s_x_map, s_y_map, holes_mask)
        elif divine_method == "gs_enhanced_divine":
            # Use full G-S coordinate information if available
            divined_image = self.gs_enhanced_divine_process(
                divined_image, original_texture, s_x_map, s_y_map, holes_mask)
        else:
            # Fallback to flow method
            divined_image = self.s_coordinate_flow_divination(
                divined_image, original_texture, s_x_map, s_y_map, holes_mask)
        
        # Count remaining holes
        remaining_holes_mask = np.all(divined_image == 0, axis=2)
        remaining_holes = np.sum(remaining_holes_mask)
        divined_holes = initial_holes - remaining_holes
        
        if self.debug_divine_process.get():
            print(f"  Successfully divined: {divined_holes}/{initial_holes} pixels")
            print(f"  Divine success rate: {(divined_holes/initial_holes)*100:.1f}%")
        
        # Final cleanup for any remaining holes
        if remaining_holes > 0:
            if self.debug_divine_process.get():
                print(f"  Applying fallback for {remaining_holes} remaining holes")
            divined_image = self.divine_fallback_cleanup(
                divined_image, original_texture, remaining_holes_mask)
        
        return divined_image

    def s_coordinate_flow_divination(self, image, original_texture, s_x, s_y, holes_mask):
        """Divine pixels by following S-coordinate flow patterns."""
        
        height, width = image.shape[:2]
        result = image.copy()
        
        # Calculate S-coordinate flow field (gradients)
        s_x_grad_x = np.gradient(s_x, axis=1)
        s_x_grad_y = np.gradient(s_x, axis=0)
        s_y_grad_x = np.gradient(s_y, axis=1)
        s_y_grad_y = np.gradient(s_y, axis=0)
        
        # Flow field vectors
        flow_x = s_x_grad_x + s_y_grad_x  # Combined X flow
        flow_y = s_x_grad_y + s_y_grad_y  # Combined Y flow
        
        # Normalize flow vectors
        flow_magnitude = np.sqrt(flow_x**2 + flow_y**2)
        flow_magnitude[flow_magnitude == 0] = 1e-6  # Avoid division by zero
        
        flow_x_norm = flow_x / flow_magnitude
        flow_y_norm = flow_y / flow_magnitude
        
        if self.debug_divine_process.get():
            print(f"  S-coordinate flow analysis:")
            print(f"    Flow magnitude range: [{flow_magnitude.min():.4f}, {flow_magnitude.max():.4f}]")
            print(f"    Flow direction variation: X[{flow_x_norm.min():.3f},{flow_x_norm.max():.3f}], Y[{flow_y_norm.min():.3f},{flow_y_norm.max():.3f}]")
        
        # For each hole, trace back along flow lines to find source pixel
        hole_coords = np.column_stack(np.where(holes_mask))
        divine_success = 0
        
        for hole_y, hole_x in hole_coords:
            # Get flow direction at hole location
            flow_dir_x = flow_x_norm[hole_y, hole_x]
            flow_dir_y = flow_y_norm[hole_y, hole_x]
            
            # Skip if flow is too weak
            if flow_magnitude[hole_y, hole_x] < self.divine_precision.get() * 0.1:
                continue
            
            # Trace back along flow line to find source
            source_pixel = self.trace_s_coordinate_flow_source(
                hole_x, hole_y, flow_dir_x, flow_dir_y, 
                s_x, s_y, original_texture, holes_mask)
            
            if source_pixel is not None:
                if len(source_pixel.shape) == 1 and len(source_pixel) == 3:
                    result[hole_y, hole_x, :] = source_pixel
                else:
                    result[hole_y, hole_x, :] = source_pixel.flatten()[:3]
                divine_success += 1
        
        if self.debug_divine_process.get():
            print(f"    Flow divination success: {divine_success}/{len(hole_coords)} holes")
        
        return result

    def trace_s_coordinate_flow_source(self, hole_x, hole_y, flow_dir_x, flow_dir_y, 
                                      s_x, s_y, original_texture, holes_mask):
        """Trace along S-coordinate flow to find the natural source pixel."""
        
        height, width = s_x.shape
        
        # Current S-coordinate values at hole
        target_s_x = s_x[hole_y, hole_x]
        target_s_y = s_y[hole_y, hole_x]
        
        # Trace backwards along flow line
        max_trace_distance = self.divine_search_radius.get()
        step_size = 0.5  # Sub-pixel tracing
        
        current_x = float(hole_x)
        current_y = float(hole_y)
        
        for step in range(int(max_trace_distance / step_size)):
            # Move backwards along flow
            current_x -= flow_dir_x * step_size
            current_y -= flow_dir_y * step_size
            
            # Check bounds
            if current_x < 0 or current_x >= width-1 or current_y < 0 or current_y >= height-1:
                break
                
            # Get integer coordinates for sampling
            sample_x = int(np.round(current_x))
            sample_y = int(np.round(current_y))
            
            # Check if this location has valid data (not a hole)
            if not holes_mask[sample_y, sample_x]:
                # Found a valid source pixel
                
                # Calculate S-coordinate similarity
                source_s_x = s_x[sample_y, sample_x]
                source_s_y = s_y[sample_y, sample_x]
                
                s_distance = np.sqrt((source_s_x - target_s_x)**2 + (source_s_y - target_s_y)**2)
                
                # If S-coordinates are similar enough, use this pixel
                similarity_threshold = self.divine_precision.get()
                if s_distance <= similarity_threshold:
                    return original_texture[sample_y, sample_x]
        
        return None

    def s_coordinate_similarity_divination(self, image, original_texture, s_x, s_y, holes_mask):
        """Divine pixels by finding S-coordinate similarity matches."""
        
        height, width = image.shape[:2]
        result = image.copy()
        
        # Build KD-tree for fast S-coordinate lookups
        valid_mask = ~holes_mask
        valid_coords = np.column_stack(np.where(valid_mask))
        
        if len(valid_coords) == 0:
            if self.debug_divine_process.get():
                print("  Warning: No valid pixels for S-coordinate similarity matching")
            return result
        
        # S-coordinates of valid pixels
        valid_s_coords = np.column_stack([
            s_x[valid_mask].ravel(),
            s_y[valid_mask].ravel()
        ])
        
        # Build KD-tree for fast nearest neighbor search in S-coordinate space
        s_coord_tree = cKDTree(valid_s_coords)
        
        if self.debug_divine_process.get():
            print(f"  S-coordinate similarity tree built with {len(valid_s_coords)} valid points")
        
        # Process each hole
        hole_coords = np.column_stack(np.where(holes_mask))
        divine_success = 0
        
        for hole_y, hole_x in hole_coords:
            # Get S-coordinates of hole
            hole_s_x = s_x[hole_y, hole_x]
            hole_s_y = s_y[hole_y, hole_x]
            hole_s_coord = np.array([hole_s_x, hole_s_y])
            
            # Find closest S-coordinate matches
            num_neighbors = min(5, len(valid_s_coords))  # Find up to 5 closest matches
            try:
                distances, indices = s_coord_tree.query(hole_s_coord, k=num_neighbors)
                
                # Check if closest match is within similarity threshold
                if distances[0] <= self.divine_precision.get():
                    # Weight pixels by S-coordinate similarity
                    weights = 1.0 / (distances + 1e-6)  # Inverse distance weighting
                    weights = weights / np.sum(weights)  # Normalize
                    
                    # Get corresponding pixel locations and values
                    divined_pixel = np.zeros(3)
                    
                    for i, idx in enumerate(indices):
                        pixel_y, pixel_x = valid_coords[idx]
                        pixel_value = original_texture[pixel_y, pixel_x]
                        divined_pixel += pixel_value * weights[i]

                    result[hole_y, hole_x, :] = divined_pixel.astype(np.uint8)
                    divine_success += 1
                    
            except Exception as e:
                if self.debug_divine_process.get():
                    print(f"    KDTree query failed for hole at ({hole_x}, {hole_y}): {e}")
                continue
        
        if self.debug_divine_process.get():
            print(f"    Similarity divination success: {divine_success}/{len(hole_coords)} holes")
        
        return result

    def s_coordinate_geometric_extrapolation(self, image, original_texture, s_x, s_y, holes_mask):
        """Divine pixels using S-coordinate geometric extrapolation."""
        
        height, width = image.shape[:2]
        result = image.copy()
        
        # Create interpolation functions for S-coordinates and texture
        valid_mask = ~holes_mask
        
        if np.sum(valid_mask) < 4:  # Need minimum points for interpolation
            if self.debug_divine_process.get():
                print("  Warning: Insufficient valid pixels for geometric extrapolation")
            return result
        
        # Get valid pixel coordinates and values
        valid_y, valid_x = np.where(valid_mask)
        valid_points = np.column_stack([valid_y, valid_x])
        
        # Valid S-coordinates
        valid_s_x = s_x[valid_mask]
        valid_s_y = s_y[valid_mask]
        
        # Valid texture values
        valid_texture = original_texture[valid_mask]
        
        if self.debug_divine_process.get():
            print(f"  Geometric extrapolation using {len(valid_points)} valid pixels")
        
        # Process holes
        hole_coords = np.column_stack(np.where(holes_mask))
        divine_success = 0
        
        for hole_y, hole_x in hole_coords:
            # Get S-coordinates at hole location
            hole_s_x = s_x[hole_y, hole_x]
            hole_s_y = s_y[hole_y, hole_x]
            
            # Find pixels with similar S-coordinates using geometric relationship
            s_distances = np.sqrt((valid_s_x - hole_s_x)**2 + (valid_s_y - hole_s_y)**2)
            
            # Use closest S-coordinate matches for interpolation
            num_closest = min(8, len(s_distances))  # Use up to 8 closest points
            closest_indices = np.argpartition(s_distances, num_closest)[:num_closest]
            
            # Check if we have close enough matches
            min_s_distance = s_distances[closest_indices[0]]
            if min_s_distance > self.divine_precision.get() * 2:  # More lenient for geometric method
                continue
            
            # Weight by S-coordinate proximity and spatial proximity
            s_weights = 1.0 / (s_distances[closest_indices] + 1e-6)
            
            # Spatial distances for additional weighting
            spatial_distances = np.sqrt((valid_y[closest_indices] - hole_y)**2 + 
                                      (valid_x[closest_indices] - hole_x)**2)
            spatial_weights = 1.0 / (spatial_distances + 1e-6)
            
            # Combined weighting (favor S-coordinate similarity over spatial proximity)
            combined_weights = s_weights * 0.7 + spatial_weights * 0.3
            combined_weights = combined_weights / np.sum(combined_weights)
            
            # Interpolate pixel value
            divined_pixel = np.average(valid_texture[closest_indices], 
                                     weights=combined_weights, axis=0)
            
            result[hole_y, hole_x, :] = divined_pixel.astype(np.uint8)
            divine_success += 1
        
        if self.debug_divine_process.get():
            print(f"    Geometric extrapolation success: {divine_success}/{len(hole_coords)} holes")
        
        return result

    def combined_divine_process(self, image, original_texture, s_x, s_y, holes_mask):
        """Multi-method divine process: try each method in sequence."""
        
        result = image.copy()
        remaining_holes = holes_mask.copy()
        
        if self.debug_divine_process.get():
            print("  🔮 Combined divine process - trying multiple methods:")
        
        # Method 1: Flow divination
        if np.any(remaining_holes):
            result = self.s_coordinate_flow_divination(result, original_texture, s_x, s_y, remaining_holes)
            remaining_holes = np.all(result == 0, axis=2)
            if self.debug_divine_process.get():
                print(f"    After flow: {np.sum(remaining_holes)} holes remaining")
        
        # Method 2: Similarity divination for remaining holes
        if np.any(remaining_holes):
            result = self.s_coordinate_similarity_divination(result, original_texture, s_x, s_y, remaining_holes)
            remaining_holes = np.all(result == 0, axis=2)
            if self.debug_divine_process.get():
                print(f"    After similarity: {np.sum(remaining_holes)} holes remaining")
        
        # Method 3: Geometric extrapolation for remaining holes
        if np.any(remaining_holes):
            result = self.s_coordinate_geometric_extrapolation(result, original_texture, s_x, s_y, remaining_holes)
            remaining_holes = np.all(result == 0, axis=2)
            if self.debug_divine_process.get():
                print(f"    After geometric: {np.sum(remaining_holes)} holes remaining")
        
        return result

    def gs_enhanced_divine_process(self, image, original_texture, s_x, s_y, holes_mask):
        """Enhanced divine process using full G-S coordinate information."""
        
        if self.gs_coordinate_map is None:
            if self.debug_divine_process.get():
                print("  No G-S coordinate maps available, falling back to combined divine")
            return self.combined_divine_process(image, original_texture, s_x, s_y, holes_mask)
        
        height, width = image.shape[:2]
        result = image.copy()
        
        # Get G-S coordinate maps
        g_theta = self.gs_coordinate_map['g_theta']
        g_phi = self.gs_coordinate_map['g_phi']
        
        # Resize if needed
        if g_theta.shape != (height, width):
            g_theta = cv2.resize(g_theta, (width, height), interpolation=cv2.INTER_LINEAR)
            g_phi = cv2.resize(g_phi, (width, height), interpolation=cv2.INTER_LINEAR)
        
        if self.debug_divine_process.get():
            print("  🌟 G-S enhanced divine process using full coordinate information")
        
        # Process holes using G-S enhanced methods
        hole_coords = np.column_stack(np.where(holes_mask))
        divine_success = 0
        
        for hole_y, hole_x in hole_coords:
            # Use comprehensive G-S coordinate synthesis
            divined_pixel = self.s_coordinate_based_pixel_synthesis(
                (hole_y, hole_x), s_x, s_y, original_texture, self.gs_coordinate_map)
            
            if divined_pixel is not None:
                if len(divined_pixel.shape) == 1 and len(divined_pixel) == 3:
                    result[hole_y, hole_x, :] = divined_pixel
                else:
                    result[hole_y, hole_x, :] = divined_pixel.flatten()[:3]
                divine_success += 1
        
        if self.debug_divine_process.get():
            print(f"    G-S enhanced divine success: {divine_success}/{len(hole_coords)} holes")
        
        # For remaining holes, use combined method
        remaining_holes = np.all(result == 0, axis=2)
        if np.any(remaining_holes):
            result = self.combined_divine_process(result, original_texture, s_x, s_y, remaining_holes)
        
        return result

    def s_coordinate_based_pixel_synthesis(self, hole_location, s_x_map, s_y_map, original_texture, 
                                         gs_coordinate_map=None):
        """
        Synthesize a pixel value based on S-coordinate geometric relationships.
        
        This is the core "divination" function that determines what pixel SHOULD exist
        at a given location based on the G-S coordinate system geometry.
        """
        hole_y, hole_x = hole_location
        height, width = s_x_map.shape
        
        # Get S-coordinates at hole location
        hole_s_x = s_x_map[hole_y, hole_x]
        hole_s_y = s_y_map[hole_y, hole_x]
        
        # Method 1: Use G-S coordinate maps for enhanced synthesis
        if gs_coordinate_map is not None:
            synthesized_pixel = self.gs_enhanced_pixel_synthesis(
                hole_s_x, hole_s_y, hole_y, hole_x, gs_coordinate_map, original_texture)
            if synthesized_pixel is not None:
                return synthesized_pixel
        
        # Method 2: S-coordinate geometric synthesis
        # Find the "natural" pixel that should exist at these S-coordinates
        
        # Calculate S-coordinate sphere position
        s_magnitude = np.sqrt(hole_s_x**2 + hole_s_y**2)
        
        if s_magnitude <= 1.0:  # Within unit sphere
            # Use sphere intersection to determine depth
            sphere_depth = np.sqrt(1.0 - s_magnitude**2)
            
            # Convert back to G-S coordinates
            g_theta_estimated = np.arctan2(hole_s_y, hole_s_x) * 180.0 / np.pi
            g_phi_estimated = np.arccos(sphere_depth) * 180.0 / np.pi
            
            # Find pixels with similar G-S coordinates
            synthesized_pixel = self.find_similar_gs_coordinate_pixel(
                g_theta_estimated, g_phi_estimated, original_texture, gs_coordinate_map)
            
            if synthesized_pixel is not None:
                return synthesized_pixel
        
        # Method 3: S-coordinate gradient-based synthesis
        # Use local S-coordinate gradients to estimate pixel value
        return self.gradient_based_pixel_synthesis(
            hole_s_x, hole_s_y, hole_y, hole_x, s_x_map, s_y_map, original_texture)

    def gs_enhanced_pixel_synthesis(self, hole_s_x, hole_s_y, hole_y, hole_x, 
                                   gs_coordinate_map, original_texture):
        """Enhanced pixel synthesis using full G-S coordinate information."""
        
        if 'g_theta' not in gs_coordinate_map or 'g_phi' not in gs_coordinate_map:
            return None
        
        height, width = gs_coordinate_map['g_theta'].shape
        
        # Get G-coordinates at hole location
        hole_g_theta = gs_coordinate_map['g_theta'][hole_y, hole_x]
        hole_g_phi = gs_coordinate_map['g_phi'][hole_y, hole_x]
        
        # Find pixels with similar G-coordinates
        g_theta_map = gs_coordinate_map['g_theta']
        g_phi_map = gs_coordinate_map['g_phi']
        
        # Calculate G-coordinate distances (handle theta wraparound)
        theta_diff = np.abs(g_theta_map - hole_g_theta)
        theta_diff = np.minimum(theta_diff, 360 - theta_diff)  # Handle wraparound
        phi_diff = np.abs(g_phi_map - hole_g_phi)
        
        # Combined G-coordinate distance
        g_distance = np.sqrt(theta_diff**2 + phi_diff**2)
        
        # Find closest matches
        similarity_threshold = self.divine_similarity_threshold.get()
        similar_mask = g_distance < similarity_threshold
        
        if np.any(similar_mask):
            # Weight by G-coordinate similarity
            weights = 1.0 / (g_distance[similar_mask] + 1e-6)
            weights = weights / np.sum(weights)
            
            # Weighted average of similar pixels
            similar_pixels = original_texture[similar_mask]
            synthesized_pixel = np.average(similar_pixels, weights=weights, axis=0)
            
            return synthesized_pixel.astype(np.uint8).flatten()
        
        return None

    def find_similar_gs_coordinate_pixel(self, target_g_theta, target_g_phi, 
                                       original_texture, gs_coordinate_map):
        """Find pixel with similar G-coordinates."""
        
        if gs_coordinate_map is None:
            return None
        
        g_theta_map = gs_coordinate_map['g_theta']
        g_phi_map = gs_coordinate_map['g_phi']
        
        # Calculate distances in G-coordinate space
        theta_diff = np.abs(g_theta_map - target_g_theta)
        theta_diff = np.minimum(theta_diff, 360 - theta_diff)  # Handle wraparound
        phi_diff = np.abs(g_phi_map - target_g_phi)
        
        g_distance = np.sqrt(theta_diff**2 + phi_diff**2)
        
        # Find closest match
        min_distance_idx = np.unravel_index(np.argmin(g_distance), g_distance.shape)
        
        if g_distance[min_distance_idx] < self.divine_similarity_threshold.get():
            return original_texture[min_distance_idx]
        
        return None

    def gradient_based_pixel_synthesis(self, hole_s_x, hole_s_y, hole_y, hole_x, 
                                     s_x_map, s_y_map, original_texture):
        """Synthesize pixel using S-coordinate gradient information."""
        
        height, width = s_x_map.shape
        
        # Calculate local S-coordinate gradients
        s_x_grad_x = np.gradient(s_x_map, axis=1)
        s_x_grad_y = np.gradient(s_x_map, axis=0)
        s_y_grad_x = np.gradient(s_y_map, axis=1)
        s_y_grad_y = np.gradient(s_y_map, axis=0)
        
        # Get gradients at hole location
        grad_sx_x = s_x_grad_x[hole_y, hole_x]
        grad_sx_y = s_x_grad_y[hole_y, hole_x]
        grad_sy_x = s_y_grad_x[hole_y, hole_x]
        grad_sy_y = s_y_grad_y[hole_y, hole_x]
        
        # Find neighboring pixels with similar gradient patterns
        search_radius = self.divine_search_radius.get()
        y_min = max(0, hole_y - search_radius)
        y_max = min(height, hole_y + search_radius + 1)
        x_min = max(0, hole_x - search_radius)
        x_max = min(width, hole_x + search_radius + 1)
        
        # Extract search region
        region_grad_sx_x = s_x_grad_x[y_min:y_max, x_min:x_max]
        region_grad_sx_y = s_x_grad_y[y_min:y_max, x_min:x_max]
        region_grad_sy_x = s_y_grad_x[y_min:y_max, x_min:x_max]
        region_grad_sy_y = s_y_grad_y[y_min:y_max, x_min:x_max]
        region_texture = original_texture[y_min:y_max, x_min:x_max]
        
        # Calculate gradient similarity
        grad_diff_sx_x = np.abs(region_grad_sx_x - grad_sx_x)
        grad_diff_sx_y = np.abs(region_grad_sx_y - grad_sx_y)
        grad_diff_sy_x = np.abs(region_grad_sy_x - grad_sy_x)
        grad_diff_sy_y = np.abs(region_grad_sy_y - grad_sy_y)
        
        # Combined gradient difference
        total_grad_diff = grad_diff_sx_x + grad_diff_sx_y + grad_diff_sy_x + grad_diff_sy_y
        
        # Find most similar gradient
        min_grad_idx = np.unravel_index(np.argmin(total_grad_diff), total_grad_diff.shape)
        
        return region_texture[min_grad_idx]

    def divine_fallback_cleanup(self, image, original_texture, holes_mask):
        """Final cleanup for any holes that couldn't be divined."""
        
        result = image.copy()
        hole_coords = np.column_stack(np.where(holes_mask))
        
        if self.debug_divine_process.get():
            print(f"  Applying fallback cleanup for {len(hole_coords)} remaining holes")
        
        # Simple nearest neighbor fallback
        for hole_y, hole_x in hole_coords:
            # Find nearest non-hole pixel
            search_radius = 1
            found = False
            
            while search_radius <= 5 and not found:
                y_min = max(0, hole_y - search_radius)
                y_max = min(image.shape[0], hole_y + search_radius + 1)
                x_min = max(0, hole_x - search_radius)
                x_max = min(image.shape[1], hole_x + search_radius + 1)
                
                for dy in range(y_min, y_max):
                    for dx in range(x_min, x_max):
                        if not np.all(result[dy, dx] == 0):
                            result[hole_y, hole_x, :] = result[dy, dx, :]
                            found = True
                            break
                    if found:
                        break
                
                if not found:
                    search_radius += 1
            
            # Ultimate fallback: use original texture
            if not found:
                result[hole_y, hole_x, :] = original_texture[hole_y, hole_x, :]
        
        return result

    # === SPHERICAL TO PLANAR CONVERSION ===
    
    def gs_spherical_to_planar_mapping(self, gs_depth_map, gs_coordinate_map=None):
        """
        Convert G-S spherical depth map to planar depth while preserving depth information.
        
        Args:
            gs_depth_map: The original G-S depth map (spherical geometry)
            gs_coordinate_map: Optional G-S coordinate maps for enhanced conversion
        
        Returns:
            planar_depth_map: Flattened depth map suitable for planar stereo
            mapping_info: Information about the conversion for debugging
        """
        height, width = gs_depth_map.shape
        center_y, center_x = height // 2, width // 2
        
        # Create coordinate grids
        y_coords, x_coords = np.mgrid[0:height, 0:width]
        
        # Calculate spherical coordinates from image position
        # Normalize to [-1, 1] range
        x_norm = (x_coords - center_x) / (width / 2)
        y_norm = (y_coords - center_y) / (height / 2)
        
        # Calculate radial distance from center (spherical coordinate)
        radius_norm = np.sqrt(x_norm**2 + y_norm**2)
        
        # Clamp radius to unit circle
        radius_norm = np.clip(radius_norm, 0, 1)
        
        # Method 1: Use G-S coordinate maps if available for precise conversion
        if gs_coordinate_map is not None and 'g_theta' in gs_coordinate_map:
            planar_depth = self.gs_coordinate_to_planar_conversion(
                gs_depth_map, gs_coordinate_map, x_norm, y_norm, radius_norm)
        else:
            # Method 2: Geometric sphere-to-plane projection
            planar_depth = self.geometric_sphere_to_plane_conversion(
                gs_depth_map, x_norm, y_norm, radius_norm)
        
        # Create mapping information for debugging
        mapping_info = {
            'original_range': [gs_depth_map.min(), gs_depth_map.max()],
            'planar_range': [planar_depth.min(), planar_depth.max()],
            'center_point': [center_x, center_y],
            'max_radius': radius_norm.max(),
            'conversion_method': 'gs_coordinate' if gs_coordinate_map else 'geometric'
        }
        
        if self.debug_divine_process.get():
            print(f"G-S to Planar Conversion:")
            print(f"  Method: {mapping_info['conversion_method']}")
            print(f"  Original depth range: [{mapping_info['original_range'][0]:.3f}, {mapping_info['original_range'][1]:.3f}]")
            print(f"  Planar depth range: [{mapping_info['planar_range'][0]:.3f}, {mapping_info['planar_range'][1]:.3f}]")
        
        return planar_depth, mapping_info

    def gs_coordinate_to_planar_conversion(self, gs_depth, gs_coords, x_norm, y_norm, radius_norm):
        """Convert using G-S coordinate information for precise mapping."""
        
        # Extract G-S coordinate maps
        g_theta = gs_coords['g_theta'] * np.pi / 180.0  # Convert to radians
        g_phi = gs_coords['g_phi'] * np.pi / 180.0
        s_x = gs_coords['s_x']
        s_y = gs_coords['s_y']
        
        # Calculate the "planar equivalent" depth
        # This preserves the relative depth relationships but maps them to planar space
        s_magnitude = np.sqrt(s_x**2 + s_y**2)
        
        # Convert from spherical depth to planar depth
        # Use the sphere equation in reverse: if depth was from sphere, what would planar depth be?
        sphere_radius = 1.0  # Normalized sphere
        
        # For points inside the sphere, calculate equivalent planar depth
        valid_mask = s_magnitude <= sphere_radius
        
        planar_depth = np.zeros_like(gs_depth)
        
        # For valid points, use the G-S depth but correct for spherical distortion
        # The idea: preserve depth ordering but remove spherical curvature
        planar_depth[valid_mask] = gs_depth[valid_mask]
        
        # Apply distortion correction based on distance from center
        # Points further from center were more "compressed" in spherical space
        distortion_correction = 1.0 + (s_magnitude * 0.3)  # 30% max correction
        distortion_correction = np.clip(distortion_correction, 1.0, 1.5)  # Reasonable limits
        
        planar_depth = planar_depth * distortion_correction
        
        # For points outside sphere (should be rare), use geometric fallback
        invalid_mask = ~valid_mask
        if np.any(invalid_mask):
            planar_depth[invalid_mask] = self.geometric_sphere_to_plane_conversion(
                gs_depth[invalid_mask], x_norm[invalid_mask], y_norm[invalid_mask], radius_norm[invalid_mask])
        
        return planar_depth

    def geometric_sphere_to_plane_conversion(self, gs_depth, x_norm, y_norm, radius_norm):
        """Geometric conversion from spherical to planar depth."""
        
        # Method: Inverse stereographic projection
        # This maps the spherical surface back to a plane while preserving depth relationships
        
        # Calculate angle from center (spherical coordinate)
        theta = radius_norm * np.pi / 2  # Map [0,1] radius to [0, π/2] angle
        
        # Stereographic projection correction
        # Points at edge of sphere (theta near π/2) were compressed
        # We need to "uncompress" them for planar viewing
        
        # Calculate correction factor
        # As theta approaches π/2, points were increasingly compressed
        correction_factor = 1.0 / np.cos(theta)
        
        # Avoid division by zero at edges
        correction_factor = np.clip(correction_factor, 1.0, 3.0)  # Max 3x correction
        
        # Apply correction to preserve depth relationships in planar space
        planar_depth = gs_depth * correction_factor
        
        # Additional radial depth correction
        # Inner points (small radius) keep original depth
        # Outer points (large radius) get depth boost to compensate for sphere curvature
        radial_boost = 1.0 + (radius_norm**2) * 0.5  # Quadratic boost based on distance from center
        
        planar_depth = planar_depth * radial_boost
        
        # Smooth the correction to avoid sharp transitions
        planar_depth = cv2.GaussianBlur(planar_depth, (3, 3), 0.5)
        
        return planar_depth

    def create_planar_stereo_displacement(self, texture, planar_depth_map, eye_sep):
        """Create stereo displacement using flattened planar depth."""
        
        height, width = texture.shape[:2]
        
        # Normalize depth to proper disparity range
        depth_normalized = (planar_depth_map - planar_depth_map.min()) / (planar_depth_map.max() - planar_depth_map.min())
        
        # Calculate disparity (standard planar stereo)
        max_disparity = self.max_disparity.get()
        disparity = (depth_normalized - 0.5) * max_disparity * eye_sep
        
        # Create coordinate grids
        y_coords, x_coords = np.mgrid[0:height, 0:width]
        
        # Calculate left and right eye positions (pure horizontal displacement)
        left_x = x_coords - disparity / 2
        right_x = x_coords + disparity / 2
        
        # Clamp to valid coordinates
        left_x = np.clip(left_x, 0, width - 1)
        right_x = np.clip(right_x, 0, width - 1)
        
        # Apply displacement
        left_image = np.zeros_like(texture)
        right_image = np.zeros_like(texture)
        
        for channel in range(3):
            left_image[:, :, channel] = self.bilinear_interpolate_safe(
                texture[:, :, channel], left_x, y_coords)
            right_image[:, :, channel] = self.bilinear_interpolate_safe(
                texture[:, :, channel], right_x, y_coords)
        
        if self.debug_divine_process.get():
            print(f"Planar stereo displacement:")
            print(f"  Disparity range: [{disparity.min():.2f}, {disparity.max():.2f}] pixels")
            print(f"  Max displacement: {max_disparity * eye_sep:.1f} pixels")
        
        return left_image, right_image

    # === MAIN DIVINE STEREO GENERATION ===
    
    def generate_divine_gs_stereo(self):
        """Generate stereo view with S-coordinate pixel divination and G-S planar conversion."""
        if self.texture_image is None or self.depth_map is None:
            messagebox.showwarning("Missing Data", "Please load texture and G-S depth map")
            return
            
        try:
            self.generate_button.config(state="disabled")
            self.progress_var.set(0)
            start_time = time.time()
            
            print("🔮 Starting Divine G-S Stereo Generation...")
            
            # Step 1: Generate G-S coordinate maps if not already done
            if self.gs_coordinate_map is None and (self.gs_divine_enabled.get() or self.gs_subpixel_enabled.get()):
                print("Generating G-S coordinate maps...")
                self.generate_gs_coordinate_maps()
            
            self.progress_var.set(20)
            self.root.update_idletasks()
            
            # Step 2: Convert G-S spherical depth to planar depth
            print("Converting G-S spherical depth to planar coordinates...")
            planar_depth_map, mapping_info = self.gs_spherical_to_planar_mapping(
                self.depth_map, self.gs_coordinate_map
            )
            
            self.progress_var.set(35)
            self.root.update_idletasks()
            
            # Step 3: Apply depth processing to planar depth
            processed_planar_depth = self.apply_depth_processing(planar_depth_map)
            
            # Step 4: Prepare texture and depth for processing
            texture_height, texture_width = self.texture_image.shape[:2]
            if processed_planar_depth.shape != (texture_height, texture_width):
                processed_planar_depth = cv2.resize(processed_planar_depth, 
                                                  (texture_width, texture_height),
                                                  interpolation=cv2.INTER_LINEAR)
            
            # Apply sub-pixel enhancement if enabled
            if self.gs_subpixel_enabled.get() and self.subpixel_enhancement_map is not None:
                enhancement_factor = int(self.gs_enhancement_factor.get())
                enhanced_texture = cv2.resize(self.texture_image, 
                                            (texture_width * enhancement_factor, texture_height * enhancement_factor),
                                            interpolation=cv2.INTER_CUBIC)
                enhanced_depth = cv2.resize(processed_planar_depth,
                                          (texture_width * enhancement_factor, texture_height * enhancement_factor),
                                          interpolation=cv2.INTER_LINEAR)
                print(f"Using sub-pixel enhancement: {enhancement_factor}x resolution")
            else:
                enhanced_texture = self.texture_image
                enhanced_depth = processed_planar_depth
            
            self.progress_var.set(50)
            self.root.update_idletasks()
            
            # Step 5: Generate stereo displacement with divine pixel process
            print("Generating divine stereo displacement...")
            displacement_method = self.displacement_method.get()
            
            if displacement_method == "gs_divine_planar":
                self.left_eye_image, self.right_eye_image = self.create_divine_planar_stereo_displacement(
                    enhanced_texture, enhanced_depth, self.eye_separation.get()
                )
            elif displacement_method == "gs_divine_spherical":
                self.left_eye_image, self.right_eye_image = self.create_divine_spherical_stereo_displacement(
                    enhanced_texture, enhanced_depth, self.eye_separation.get()
                )
            else:
                # Fallback to standard enhanced displacement
                self.left_eye_image, self.right_eye_image = self.create_planar_stereo_displacement(
                    enhanced_texture, enhanced_depth, self.eye_separation.get()
                )
            
            self.progress_var.set(75)
            self.root.update_idletasks()
            
            # Step 6: Downsample if sub-pixel enhancement was used
            if self.gs_subpixel_enabled.get() and self.subpixel_enhancement_map is not None:
                original_height, original_width = self.texture_image.shape[:2]
                self.left_eye_image = cv2.resize(self.left_eye_image, (original_width, original_height),
                                               interpolation=cv2.INTER_AREA)
                self.right_eye_image = cv2.resize(self.right_eye_image, (original_width, original_height),
                                                interpolation=cv2.INTER_AREA)
            
            # Step 7: Create combined view
            self.stereo_combined = self.create_combined_stereo_view(
                self.left_eye_image, self.right_eye_image
            )
            
            processing_time = time.time() - start_time
            
            # Update UI
            self.process_status_label.config(
                text=f"✅ Divine G-S stereo generated ({processing_time:.2f}s)", 
                foreground="green"
            )
            self.export_button.config(state="normal")
            self.export_divine_maps_button.config(state="normal")
            
            # Store processing info
            self.last_processing_info = {
                'processing_time': processing_time,
                'mapping_info': mapping_info,
                'divine_enabled': self.gs_divine_enabled.get(),
                'divine_method': self.divine_method.get(),
                'subpixel_enabled': self.gs_subpixel_enabled.get(),
                'enhancement_factor': self.gs_enhancement_factor.get() if self.gs_subpixel_enabled.get() else 1.0
            }
            
            self.update_statistics()
            self.visualize_stereo_result()
            
            print(f"🔮 Divine G-S stereo generation completed in {processing_time:.2f} seconds")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate divine G-S stereo: {str(e)}")
            print(f"Error details: {str(e)}")
        finally:
            self.generate_button.config(state="normal")
            self.progress_var.set(100)

    def create_divine_planar_stereo_displacement(self, texture, planar_depth_map, eye_sep):
        """Create stereo displacement using divine pixel process with planar depth."""
        
        height, width = texture.shape[:2]
        
        # Standard planar stereo calculation
        depth_normalized = (planar_depth_map - planar_depth_map.min()) / (planar_depth_map.max() - planar_depth_map.min())
        disparity = (depth_normalized - 0.5) * self.max_disparity.get() * eye_sep
        
        # Create coordinate grids
        y_coords, x_coords = np.mgrid[0:height, 0:width]
        
        # Calculate left and right eye positions (pure horizontal displacement)
        left_x = x_coords - disparity / 2
        right_x = x_coords + disparity / 2
        
        # Apply displacement with bounds checking
        left_image = np.zeros_like(texture)
        right_image = np.zeros_like(texture)
        
        for channel in range(3):
            left_image[:, :, channel] = self.bilinear_interpolate_safe(
                texture[:, :, channel], left_x, y_coords)
            right_image[:, :, channel] = self.bilinear_interpolate_safe(
                texture[:, :, channel], right_x, y_coords)
        
        # Apply divine pixel process if enabled
        if self.gs_divine_enabled.get():
            left_holes = np.all(left_image == 0, axis=2)
            right_holes = np.all(right_image == 0, axis=2)
            
            initial_left_holes = np.sum(left_holes)
            initial_right_holes = np.sum(right_holes)
            
            print(f"🔮 Divine process: {initial_left_holes} left holes, {initial_right_holes} right holes")
            
            # Divine pixels using S-coordinate information
            if self.gs_coordinate_map is not None and 's_x' in self.gs_coordinate_map:
                s_x_map = self.gs_coordinate_map['s_x']
                s_y_map = self.gs_coordinate_map['s_y']
                
                # Resize S-coordinate maps if needed
                if s_x_map.shape != (height, width):
                    s_x_map = cv2.resize(s_x_map, (width, height), interpolation=cv2.INTER_LINEAR)
                    s_y_map = cv2.resize(s_y_map, (width, height), interpolation=cv2.INTER_LINEAR)
                
                # Divine pixels for left eye
                if np.any(left_holes):
                    left_image = self.divine_pixels_from_s_coordinates(
                        left_image, texture, s_x_map, s_y_map, left_holes)
                
                # Divine pixels for right eye  
                if np.any(right_holes):
                    right_image = self.divine_pixels_from_s_coordinates(
                        right_image, texture, s_x_map, s_y_map, right_holes)
                
                # Report divine success
                final_left_holes = np.sum(np.all(left_image == 0, axis=2))
                final_right_holes = np.sum(np.all(right_image == 0, axis=2))
                
                left_divined = initial_left_holes - final_left_holes
                right_divined = initial_right_holes - final_right_holes
                
                print(f"🔮 Divine success: Left {left_divined}/{initial_left_holes}, Right {right_divined}/{initial_right_holes}")
            else:
                print("⚠️ Warning: No S-coordinate maps available for divine process")
                # Fallback to simple hole filling
                left_image[left_holes, :] = texture[left_holes, :]
                right_image[right_holes, :] = texture[right_holes, :]
        
        return left_image, right_image

    def create_divine_spherical_stereo_displacement(self, texture, depth_map, eye_sep):
        """Create stereo displacement preserving spherical geometry with divine pixels."""
        
        height, width = texture.shape[:2]
        center_y, center_x = height // 2, width // 2
        
        # Calculate spherical displacement (preserves G-S geometry)
        y_coords, x_coords = np.mgrid[0:height, 0:width]
        
        # Distance from center (spherical coordinate)
        radius_from_center = np.sqrt((y_coords - center_y)**2 + (x_coords - center_x)**2)
        max_radius = min(center_y, center_x)
        
        # Spherical correction factor
        sphere_factor = np.cos(radius_from_center / max_radius * np.pi / 2)
        sphere_factor = np.clip(sphere_factor, 0.1, 1.0)
        
        # Apply depth with spherical compensation
        disparity = (depth_map - 0.5) * self.max_disparity.get() * eye_sep * sphere_factor
        
        # Standard displacement
        left_x = x_coords - disparity / 2
        right_x = x_coords + disparity / 2
        
        # Apply displacement
        left_image = np.zeros_like(texture)
        right_image = np.zeros_like(texture)
        
        for channel in range(3):
            left_image[:, :, channel] = self.bilinear_interpolate_safe(
                texture[:, :, channel], left_x, y_coords)
            right_image[:, :, channel] = self.bilinear_interpolate_safe(
                texture[:, :, channel], right_x, y_coords)
        
        # Apply divine process as in planar version
        if self.gs_divine_enabled.get() and self.gs_coordinate_map is not None:
            left_holes = np.all(left_image == 0, axis=2)
            right_holes = np.all(right_image == 0, axis=2)
            
            if np.any(left_holes) or np.any(right_holes):
                s_x_map = self.gs_coordinate_map['s_x']
                s_y_map = self.gs_coordinate_map['s_y']
                
                if s_x_map.shape != (height, width):
                    s_x_map = cv2.resize(s_x_map, (width, height), interpolation=cv2.INTER_LINEAR)
                    s_y_map = cv2.resize(s_y_map, (width, height), interpolation=cv2.INTER_LINEAR)
                
                if np.any(left_holes):
                    left_image = self.divine_pixels_from_s_coordinates(
                        left_image, texture, s_x_map, s_y_map, left_holes)
                
                if np.any(right_holes):
                    right_image = self.divine_pixels_from_s_coordinates(
                        right_image, texture, s_x_map, s_y_map, right_holes)
        
        return left_image, right_image

    def test_divine_process(self):
        """Test the divine process on current data."""
        if self.texture_image is None or self.depth_map is None:
            messagebox.showwarning("Missing Data", "Please load texture and G-S depth map first")
            return
        
        if self.gs_coordinate_map is None:
            messagebox.showwarning("Missing Coordinates", "Generate G-S coordinate maps first")
            return
        
        try:
            print("🔮 Testing Divine Process...")
            
            # Create a small test area with artificial holes
            test_image = self.texture_image.copy()
            height, width = test_image.shape[:2]
            
            # Create test holes (10x10 squares)
            test_holes = np.zeros((height, width), dtype=bool)
            for i in range(3):
                y_start = height // 4 + i * height // 6
                x_start = width // 4 + i * width // 6
                test_holes[y_start:y_start+10, x_start:x_start+10] = True
                test_image[y_start:y_start+10, x_start:x_start+10] = 0
            
            s_x_map = self.gs_coordinate_map['s_x']
            s_y_map = self.gs_coordinate_map['s_y']
            
            # Resize if needed
            if s_x_map.shape != (height, width):
                s_x_map = cv2.resize(s_x_map, (width, height), interpolation=cv2.INTER_LINEAR)
                s_y_map = cv2.resize(s_y_map, (width, height), interpolation=cv2.INTER_LINEAR)
            
            # Apply divine process
            divined_image = self.divine_pixels_from_s_coordinates(
                test_image, self.texture_image, s_x_map, s_y_map, test_holes)
            
            # Visualize results
            self.visualize_divine_test_results(self.texture_image, test_image, divined_image, test_holes)
            
            print("🔮 Divine process test completed")
            
        except Exception as e:
            messagebox.showerror("Test Error", f"Divine process test failed: {str(e)}")

    def bilinear_interpolate_safe(self, image, x_coords, y_coords):
        """Enhanced safe bilinear interpolation with better bounds checking."""
        height, width = image.shape
        
        # More aggressive coordinate clamping
        x_coords = np.clip(x_coords, 0, width - 1)
        y_coords = np.clip(y_coords, 0, height - 1)
        
        # Get integer coordinates
        x0 = np.floor(x_coords).astype(int)
        x1 = np.clip(x0 + 1, 0, width - 1)
        y0 = np.floor(y_coords).astype(int)
        y1 = np.clip(y0 + 1, 0, height - 1)
        
        # Ensure all indices are valid
        x0 = np.clip(x0, 0, width - 1)
        x1 = np.clip(x1, 0, width - 1)
        y0 = np.clip(y0, 0, height - 1)
        y1 = np.clip(y1, 0, height - 1)
        
        # Get fractional parts
        wx = x_coords - x0
        wy = y_coords - y0
        
        # Clamp weights to [0,1] range
        wx = np.clip(wx, 0, 1)
        wy = np.clip(wy, 0, 1)
        
        try:
            # Bilinear interpolation with error handling
            result = (image[y0, x0] * (1 - wx) * (1 - wy) +
                     image[y0, x1] * wx * (1 - wy) +
                     image[y1, x0] * (1 - wx) * wy +
                     image[y1, x1] * wx * wy)
            
            # Check for NaN or inf values
            if np.any(np.isnan(result)) or np.any(np.isinf(result)):
                if self.debug_divine_process.get():
                    print("Warning: Invalid values in interpolation result, using nearest neighbor fallback")
                result = image[y0, x0]  # Fallback to nearest neighbor
            
            return result
        
        except IndexError as e:
            if self.debug_divine_process.get():
                print(f"Index error in bilinear interpolation: {e}")
            # Fallback to nearest neighbor
            return image[y0, x0]

    def create_combined_stereo_view(self, left_image, right_image):
        """Create combined stereo view."""
        mode = self.stereo_mode.get()
        
        if mode == "side_by_side":
            combined = np.hstack([left_image, right_image])
        elif mode == "over_under":
            combined = np.vstack([left_image, right_image])
        elif mode == "cross_eyed":
            combined = np.hstack([right_image, left_image])
        elif mode == "anaglyph_red_cyan":
            left_gray = cv2.cvtColor(left_image, cv2.COLOR_RGB2GRAY)
            right_gray = cv2.cvtColor(right_image, cv2.COLOR_RGB2GRAY)
            combined = np.zeros_like(left_image)
            combined[:, :, 0] = left_gray
            combined[:, :, 1] = right_gray
            combined[:, :, 2] = right_gray
        elif mode == "left_only":
            combined = left_image
        elif mode == "right_only":
            combined = right_image
        else:
            combined = np.hstack([left_image, right_image])
        
        # Add overlays if enabled
        if self.show_gs_overlay.get() and self.gs_coordinate_map is not None:
            combined = self.add_gs_coordinate_overlay(combined, mode)
        
        if self.show_divine_overlay.get() and self.gs_divine_enabled.get():
            combined = self.add_divine_process_overlay(combined, mode)
        
        if self.show_guides.get() and mode in ["side_by_side", "over_under", "cross_eyed"]:
            combined = self.add_alignment_guides_improved(combined, mode)
            
        return combined

    def add_gs_coordinate_overlay(self, image, mode):
        """Add G-S coordinate information overlay."""
        overlay_image = image.copy()
        
        if mode in ["side_by_side", "cross_eyed"]:
            height, width = overlay_image.shape[:2]
            overlay_width = width // 2
            
            if 'g_theta' in self.gs_coordinate_map:
                g_theta_sample = self.gs_coordinate_map['g_theta'][height//4, overlay_width//4]
                g_phi_sample = self.gs_coordinate_map['g_phi'][height//4, overlay_width//4]
                
                cv2.putText(overlay_image, f"G_theta: {g_theta_sample:.1f}°", 
                           (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 1)
                cv2.putText(overlay_image, f"G_phi: {g_phi_sample:.1f}°", 
                           (10, 50), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 1)
                cv2.putText(overlay_image, "G-S Enhanced", 
                           (10, 70), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 255), 1)
        
        return overlay_image

    def add_divine_process_overlay(self, image, mode):
        """Add divine process information overlay."""
        overlay_image = image.copy()
        
        if mode in ["side_by_side", "cross_eyed"]:
            height, width = overlay_image.shape[:2]
            
            # Add divine process info
            cv2.putText(overlay_image, f"Divine: {self.divine_method.get()}", 
                       (10, height - 70), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 0, 255), 1)
            cv2.putText(overlay_image, f"Precision: {self.divine_precision.get():.3f}", 
                       (10, height - 50), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 0, 255), 1)
            cv2.putText(overlay_image, f"Search: {self.divine_search_radius.get()}px", 
                       (10, height - 30), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 0, 255), 1)
            cv2.putText(overlay_image, "🔮 DIVINE PIXELS", 
                       (10, height - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 0, 255), 1)
        
        return overlay_image

    def add_alignment_guides_improved(self, image, mode):
        """Add subtle alignment guides."""
        guide_image = image.copy()
        height, width = guide_image.shape[:2]
        guide_color = (200, 200, 200)
        
        if mode == "side_by_side" or mode == "cross_eyed":
            center_x = width // 2
            cv2.line(guide_image, (center_x, height//4), (center_x, 3*height//4), guide_color, 1)
            
            marker_size = 20
            cv2.line(guide_image, (marker_size, marker_size), (0, marker_size), guide_color, 1)
            cv2.line(guide_image, (marker_size, marker_size), (marker_size, 0), guide_color, 1)
            cv2.line(guide_image, (width-marker_size, marker_size), (width, marker_size), guide_color, 1)
            cv2.line(guide_image, (width-marker_size, marker_size), (width-marker_size, 0), guide_color, 1)
            
        elif mode == "over_under":
            center_y = height // 2
            cv2.line(guide_image, (width//4, center_y), (3*width//4, center_y), guide_color, 1)
        
        return guide_image

    # === VISUALIZATION METHODS ===
    
    def visualize_gs_coordinate_maps(self):
        """Visualize the generated G-S coordinate maps."""
        if self.gs_coordinate_map is None:
            return
            
        self.fig.clear()
        
        # Create subplot layout
        gs = self.fig.add_gridspec(2, 3, hspace=0.4, wspace=0.3)
        
        # Original depth map
        ax1 = self.fig.add_subplot(gs[0, 0])
        im1 = ax1.imshow(self.depth_map, cmap='viridis')
        ax1.set_title("Original G-S Depth Map")
        ax1.axis('off')
        self.fig.colorbar(im1, ax=ax1, fraction=0.046)
        
        # G_theta map
        ax2 = self.fig.add_subplot(gs[0, 1])
        im2 = ax2.imshow(self.gs_coordinate_map['g_theta'], cmap='hsv')
        ax2.set_title("G_theta (degrees)")
        ax2.axis('off')
        self.fig.colorbar(im2, ax=ax2, fraction=0.046)
        
        # G_phi map
        ax3 = self.fig.add_subplot(gs[0, 2])
        im3 = ax3.imshow(self.gs_coordinate_map['g_phi'], cmap='plasma')
        ax3.set_title("G_phi (degrees)")
        ax3.axis('off')
        self.fig.colorbar(im3, ax=ax3, fraction=0.046)
        
        # S_x map
        ax4 = self.fig.add_subplot(gs[1, 0])
        im4 = ax4.imshow(self.gs_coordinate_map['s_x'], cmap='coolwarm')
        ax4.set_title("S_x coordinates")
        ax4.axis('off')
        self.fig.colorbar(im4, ax=ax4, fraction=0.046)
        
        # S_y map
        ax5 = self.fig.add_subplot(gs[1, 1])
        im5 = ax5.imshow(self.gs_coordinate_map['s_y'], cmap='coolwarm')
        ax5.set_title("S_y coordinates")
        ax5.axis('off')
        self.fig.colorbar(im5, ax=ax5, fraction=0.046)
        
        # Sub-pixel enhancement (if available)
        if self.subpixel_enhancement_map is not None:
            ax6 = self.fig.add_subplot(gs[1, 2])
            enhancement_factor = int(self.gs_enhancement_factor.get())
            display_enhancement = self.subpixel_enhancement_map[::enhancement_factor, ::enhancement_factor]
            im6 = ax6.imshow(display_enhancement, cmap='viridis')
            ax6.set_title(f"Sub-pixel Enhanced ({enhancement_factor}x)")
            ax6.axis('off')
            self.fig.colorbar(im6, ax=ax6, fraction=0.046)
        else:
            ax6 = self.fig.add_subplot(gs[1, 2])
            ax6.text(0.5, 0.5, "Enable G-S\nSub-pixel\nEnhancement", 
                    ha='center', va='center', transform=ax6.transAxes)
            ax6.set_title("Sub-pixel Enhancement")
            ax6.axis('off')
        
        self.fig.suptitle("G-S Coordinate Maps and Sub-pixel Enhancement", fontsize=14, weight='bold')
        self.canvas.draw()

    def visualize_stereo_result(self):
        """Visualize the divine G-S enhanced stereo result."""
        if self.stereo_combined is None:
            return
            
        self.fig.clear()
        
        mode = self.stereo_mode.get()
        
        if mode in ["side_by_side", "cross_eyed", "over_under"]:
            gs = self.fig.add_gridspec(2, 2, hspace=0.3, wspace=0.2)
            
            ax1 = self.fig.add_subplot(gs[0, 0])
            ax1.imshow(self.left_eye_image)
            ax1.set_title("Divine Enhanced Left Eye")
            ax1.axis('off')
            
            ax2 = self.fig.add_subplot(gs[0, 1])
            ax2.imshow(self.right_eye_image)
            ax2.set_title("Divine Enhanced Right Eye")
            ax2.axis('off')
            
            ax3 = self.fig.add_subplot(gs[1, :])
            ax3.imshow(self.stereo_combined)
            ax3.set_title(f"🔮 Divine G-S Stereo ({mode.replace('_', ' ').title()})")
            ax3.axis('off')
            
        else:
            ax = self.fig.add_subplot(111)
            ax.imshow(self.stereo_combined)
            ax.set_title(f"🔮 Divine G-S Stereo ({mode.replace('_', ' ').title()})")
            ax.axis('off')
        
        enhancement_text = ""
        if self.gs_divine_enabled.get():
            enhancement_text += f" • 🔮 Divine: {self.divine_method.get()}"
        if self.gs_subpixel_enabled.get():
            enhancement_text += f" • {self.gs_enhancement_factor.get():.1f}x Sub-pixel"
        if self.gs_coordinate_map is not None:
            enhancement_text += " • G-S Coordinates"
        if self.gs_alignment_data is not None:
            enhancement_text += f" • {len(self.gs_alignment_data)} Alignments"
        
        self.fig.suptitle(f"Divine G-S Enhanced Stereo Result{enhancement_text}", 
                         fontsize=14, weight='bold')
        self.canvas.draw()

    def visualize_divine_test_results(self, original, with_holes, divined, holes_mask):
        """Visualize the results of divine process testing."""
        
        self.fig.clear()
        
        # Create comparison layout
        gs = self.fig.add_gridspec(2, 2, hspace=0.3, wspace=0.2)
        
        # Original image
        ax1 = self.fig.add_subplot(gs[0, 0])
        ax1.imshow(original)
        ax1.set_title("Original Texture")
        ax1.axis('off')
        
        # Image with artificial holes
        ax2 = self.fig.add_subplot(gs[0, 1])
        ax2.imshow(with_holes)
        ax2.set_title("With Test Holes")
        ax2.axis('off')
        
        # Divined result
        ax3 = self.fig.add_subplot(gs[1, 0])
        ax3.imshow(divined)
        ax3.set_title("After Divine Process")
        ax3.axis('off')
        
        # Difference visualization
        ax4 = self.fig.add_subplot(gs[1, 1])
        difference = np.abs(original.astype(float) - divined.astype(float))
        difference_gray = np.mean(difference, axis=2)
        im4 = ax4.imshow(difference_gray, cmap='hot', vmin=0, vmax=50)
        ax4.set_title("Difference (Divine vs Original)")
        ax4.axis('off')
        self.fig.colorbar(im4, ax=ax4, fraction=0.046, label='Pixel Difference')
        
        # Add divine method info
        method_text = f"Divine Method: {self.divine_method.get()}\n"
        method_text += f"Precision: {self.divine_precision.get():.3f}\n"
        method_text += f"Search Radius: {self.divine_search_radius.get()}\n"
        method_text += f"Similarity Threshold: {self.divine_similarity_threshold.get():.1f}°"
        
        ax3.text(0.02, 0.98, method_text, transform=ax3.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.8))
        
        self.fig.suptitle("🔮 Divine Process Test Results", fontsize=14, weight='bold')
        self.canvas.draw()

    def update_statistics(self):
        """Update divine enhancement statistics."""
        if self.stereo_combined is None:
            return
            
        stats_text = "🔮 DIVINE G-S STATISTICS\n"
        stats_text += "=" * 35 + "\n"
        
        if self.texture_image is not None:
            stats_text += f"Texture: {self.texture_image.shape}\n"
        
        if self.depth_map is not None:
            stats_text += f"G-S Depth: {self.depth_map.shape}\n"
            stats_text += f"Range: [{self.depth_map.min():.3f}, {self.depth_map.max():.3f}]\n"
        
        if self.gs_alignment_data is not None:
            stats_text += f"Alignments: {len(self.gs_alignment_data)}\n"
        
        if self.gs_divine_enabled.get():
            stats_text += f"🔮 Divine: {self.divine_method.get()}\n"
            stats_text += f"Precision: {self.divine_precision.get():.3f}\n"
            stats_text += f"Search: {self.divine_search_radius.get()}px\n"
            stats_text += f"Threshold: {self.divine_similarity_threshold.get():.1f}°\n"
        
        if self.gs_subpixel_enabled.get():
            stats_text += f"Sub-pixel: {self.gs_enhancement_factor.get():.1f}x\n"
            stats_text += f"Interpolation: {self.gs_interpolation_method.get()}\n"
        
        if self.gs_coordinate_map is not None:
            stats_text += "G-S Coordinates: ACTIVE\n"
        
        stats_text += f"Displacement: {self.displacement_method.get()}\n"
        stats_text += f"Hole Replace: {self.hole_fill_method.get()}\n"
        stats_text += f"Eye Sep: {self.eye_separation.get():.1f}px\n"
        stats_text += f"Max Disp: {self.max_disparity.get():.1f}px\n"
        stats_text += f"Mode: {self.stereo_mode.get()}\n"
        stats_text += f"Output: {self.stereo_combined.shape}\n"
        stats_text += f"Generated: {datetime.now().strftime('%H:%M:%S')}"
        
        self.stats_label.config(text=stats_text)

    def show_loaded_data(self):
        """Show currently loaded G-S data."""
        if self.texture_image is None and self.depth_map is None:
            return
            
        self.fig.clear()
        
        if self.texture_image is not None and self.depth_map is not None:
            ax1 = self.fig.add_subplot(121)
            ax1.imshow(self.texture_image)
            ax1.set_title(f"Texture Image ({self.texture_image.shape[1]}x{self.texture_image.shape[0]})")
            ax1.axis('off')
            
            ax2 = self.fig.add_subplot(122)
            im = ax2.imshow(self.depth_map, cmap='viridis')
            ax2.set_title(f"G-S Depth Map ({self.depth_map.shape[1]}x{self.depth_map.shape[0]})")
            ax2.axis('off')
            self.fig.colorbar(im, ax=ax2, label='G-S Depth Value')
            
            if self.gs_alignment_data:
                ax2.text(0.02, 0.98, f"Alignments: {len(self.gs_alignment_data)}", 
                        transform=ax2.transAxes, verticalalignment='top',
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
            
        elif self.texture_image is not None:
            ax = self.fig.add_subplot(111)
            ax.imshow(self.texture_image)
            ax.set_title(f"Texture Image ({self.texture_image.shape[1]}x{self.texture_image.shape[0]})")
            ax.axis('off')
            
        elif self.depth_map is not None:
            ax = self.fig.add_subplot(111)
            im = ax.imshow(self.depth_map, cmap='viridis')
            ax.set_title(f"G-S Depth Map ({self.depth_map.shape[1]}x{self.depth_map.shape[0]})")
            ax.axis('off')
            self.fig.colorbar(im, ax=ax, label='G-S Depth Value')
        
        self.canvas.draw()

    # === EXPORT METHODS ===
    
    def export_divine_stereo(self):
        """Export divine enhanced stereo images."""
        if self.left_eye_image is None or self.right_eye_image is None:
            messagebox.showwarning("No Data", "Generate divine stereo first")
            return
        
        directory = filedialog.askdirectory(title="Select directory to save divine stereo")
        if not directory:
            return
        
        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            enhancement_suffix = "_divine"
            
            if self.gs_divine_enabled.get():
                enhancement_suffix += f"_{self.divine_method.get()}"
            if self.gs_subpixel_enabled.get():
                enhancement_suffix += f"_subpixel{self.gs_enhancement_factor.get():.0f}x"
            if self.gs_coordinate_map is not None:
                enhancement_suffix += "_gscoords"
            
            # Save individual eyes
            left_path = os.path.join(directory, f"divine_left{enhancement_suffix}_{timestamp}.png")
            left_pil = Image.fromarray(self.left_eye_image.astype(np.uint8))
            left_pil.save(left_path)
            
            right_path = os.path.join(directory, f"divine_right{enhancement_suffix}_{timestamp}.png")
            right_pil = Image.fromarray(self.right_eye_image.astype(np.uint8))
            right_pil.save(right_path)
            
            # Save combined view
            combined_path = os.path.join(directory, f"divine_combined{enhancement_suffix}_{timestamp}.png")
            combined_pil = Image.fromarray(self.stereo_combined.astype(np.uint8))
            combined_pil.save(combined_path)
            
            # Save metadata
            metadata_path = os.path.join(directory, f"divine_metadata{enhancement_suffix}_{timestamp}.txt")
            with open(metadata_path, 'w') as f:
                f.write(f"🔮 Divine G-S Stereo Export\n")
                f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Author: Angledcrystals\n")
                f.write(f"Version: G-S Divine Stereo Viewer v1.3\n\n")
                
                f.write(f"🔮 DIVINE PROCESS SETTINGS:\n")
                f.write(f"Divine Enabled: {'YES' if self.gs_divine_enabled.get() else 'NO'}\n")
                if self.gs_divine_enabled.get():
                    f.write(f"Divine Method: {self.divine_method.get()}\n")
                    f.write(f"Divine Precision: {self.divine_precision.get():.3f}\n")
                    f.write(f"Search Radius: {self.divine_search_radius.get()}px\n")
                    f.write(f"Similarity Threshold: {self.divine_similarity_threshold.get():.1f}°\n")
                
                f.write(f"\nG-S ENHANCEMENT SETTINGS:\n")
                f.write(f"Sub-pixel Enhancement: {'ON' if self.gs_subpixel_enabled.get() else 'OFF'}\n")
                if self.gs_subpixel_enabled.get():
                    f.write(f"Enhancement Factor: {self.gs_enhancement_factor.get():.1f}x\n")
                    f.write(f"Interpolation Method: {self.gs_interpolation_method.get()}\n")
                
                f.write(f"G-S Coordinate Maps: {'ACTIVE' if self.gs_coordinate_map else 'NONE'}\n")
                f.write(f"G-S Alignment Data: {len(self.gs_alignment_data) if self.gs_alignment_data else 'NONE'}\n")
                
                f.write(f"\nSTEREO PARAMETERS:\n")
                f.write(f"Mode: {self.stereo_mode.get()}\n")
                f.write(f"Displacement Method: {self.displacement_method.get()}\n")
                f.write(f"Hole Replacement: {self.hole_fill_method.get()}\n")
                f.write(f"Eye Separation: {self.eye_separation.get():.1f}px\n")
                f.write(f"Max Disparity: {self.max_disparity.get():.1f}px\n")
                
                f.write(f"\nOUTPUT FILES:\n")
                f.write(f"Left Eye: {os.path.basename(left_path)}\n")
                f.write(f"Right Eye: {os.path.basename(right_path)}\n")
                f.write(f"Combined: {os.path.basename(combined_path)}\n")
                f.write(f"Final Size: {self.stereo_combined.shape}\n")
            
            messagebox.showinfo("Export Success", 
                              f"🔮 Divine Stereo exported:\n" + 
                              f"• {os.path.basename(left_path)}\n" +
                              f"• {os.path.basename(right_path)}\n" +
                              f"• {os.path.basename(combined_path)}\n" +
                              f"• {os.path.basename(metadata_path)}")
            
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export divine stereo: {str(e)}")

    def export_divine_maps(self):
        """Export divine process maps and coordinate data."""
        if self.gs_coordinate_map is None:
            messagebox.showwarning("No Data", "Generate G-S coordinate maps first")
            return
        
        directory = filedialog.askdirectory(title="Select directory to save divine maps")
        if not directory:
            return
        
        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Export coordinate maps
            for coord_name, coord_map in self.gs_coordinate_map.items():
                map_path = os.path.join(directory, f"divine_map_{coord_name}_{timestamp}.png")
                
                # Normalize to 0-255 range
                normalized_map = ((coord_map - coord_map.min()) / (coord_map.max() - coord_map.min()) * 255).astype(np.uint8)
                
                map_pil = Image.fromarray(normalized_map, mode='L')
                map_pil.save(map_path)
            
            # Export raw coordinate data
            coords_path = os.path.join(directory, f"divine_coordinates_{timestamp}.npz")
            np.savez_compressed(coords_path, **self.gs_coordinate_map)
            
            messagebox.showinfo("Export Success", f"🔮 Divine maps exported to:\n{directory}")
            
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export divine maps: {str(e)}")

    # === EVENT HANDLERS ===
    
    def on_parameter_change(self, *args):
        """Handle parameter changes."""
        if (self.auto_update_var.get() and self.texture_image is not None and 
            self.depth_map is not None and self.left_eye_image is not None):
            self.root.after(300, self.generate_divine_gs_stereo)

    def update_data_status(self):
        """Update data loading status."""
        status_parts = []
        
        if self.texture_image is not None:
            status_parts.append(f"✅ Texture: {self.texture_image.shape}")
        else:
            status_parts.append("❌ No texture")
            
        if self.depth_map is not None:
            status_parts.append(f"✅ Depth: {self.depth_map.shape}")
        else:
            status_parts.append("❌ No depth map")
            
        if self.gs_alignment_data is not None:
            status_parts.append(f"✅ Alignments: {len(self.gs_alignment_data)}")
        else:
            status_parts.append("❌ No alignments")
        
        status_text = " | ".join(status_parts)
        
        if self.texture_image is not None and self.depth_map is not None:
            self.data_status_label.config(text=status_text, foreground="green")
            self.generate_button.config(state="normal")
            self.process_status_label.config(text="Ready for divine G-S stereo", foreground="blue")
        else:
            self.data_status_label.config(text=status_text, foreground="orange")

# === MAIN APPLICATION ENTRY POINT ===

def main():
    """Main application entry point."""
    print("🔮 Starting G-S Divine Stereo Viewer v1.3")
    print("S-Coordinate Pixel Divination System")
    print("Author: Angledcrystals")
    print("Date: 2025-06-09 07:28:09 UTC")
    print("=" * 50)
    
    try:
        # Create and run the application
        root = tk.Tk()
        app = GSDivineStereoViewer(root)
        
        # Set window icon if available
        try:
            root.iconname("G-S Divine Stereo")
        except:
            pass
        
        # Center window on screen
        root.update_idletasks()
        width = root.winfo_width()
        height = root.winfo_height()
        x = (root.winfo_screenwidth() // 2) - (width // 2)
        y = (root.winfo_screenheight() // 2) - (height // 2)
        root.geometry(f"{width}x{height}+{x}+{y}")
        
        print("🔮 G-S Divine Stereo Viewer launched successfully!")
        print("Ready for S-coordinate pixel divination...")
        
        # Start the main loop
        root.mainloop()
        
    except Exception as e:
        print(f"❌ Failed to start G-S Divine Stereo Viewer: {e}")
        messagebox.showerror("Startup Error", f"Failed to start application: {str(e)}")
    
    print("🔮 G-S Divine Stereo Viewer closed")

if __name__ == "__main__":
    main()
